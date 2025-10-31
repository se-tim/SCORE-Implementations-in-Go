package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"runtime/debug"
	"strings"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func main() {
	runtime.GOMAXPROCS(1) // Only one thread for reproducibility

	// Parameters

	logN := 12
	numLevelsAfterBoot := 1
	logSecretWeight := 6
	
	var numBootRuns int
	fmt.Print("Enter the number of runs: ")
	fmt.Scanf("%d", &numBootRuns)

	var logNumSlots int
	fmt.Printf("Enter the number of slots (between 2**1 and 2**%d): 2**", logN-logSecretWeight-2)
	fmt.Scanf("%d", &logNumSlots)
	if logNumSlots > logN-logSecretWeight-2 {logNumSlots = logN - logSecretWeight - 2}
	if logNumSlots < 1 {logNumSlots = 1}

	N := 1 << logN
	n := 1 << logNumSlots
	secretWeight := 1 << logSecretWeight
	B := N / secretWeight
	Bover4n := B / (4 * n)
	Bover2n := 2 * Bover4n

	logDefaultScale := 40
	logBootScale := 55
	defaultScale := rlwe.NewScale(math.Exp2(float64(logDefaultScale)))
	bootScale := rlwe.NewScale(math.Exp2(float64(logBootScale)))

	logQBase := []int{53}
	logQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range logQiAfterBoot {logQiAfterBoot[i] = logDefaultScale}
	logQiSlotsToCoeffs := []int{39}
	logQiRemaining := make([]int, logSecretWeight+1) // For product tree and cs products
	for i := range logQiRemaining {logQiRemaining[i] = logBootScale}

	logQ := append(logQBase, logQiAfterBoot...)
	logQ = append(logQ, logQiSlotsToCoeffs...)
	logQ = append(logQ, logQiRemaining...)

	logP := []int{61}

	logPQ := 0
	for _, Q := range logQ {logPQ += Q}
	for _, P := range logP {logPQ += P}

	println()
	fmt.Printf("logN = %d, logNumSlots = %d, logPQ = %d, logQBeforeBoot = %d, logQAfterBoot = %d\n\n",
		logN,
		logNumSlots,
		logPQ,
		logQBase[0],
		logQBase[0] + logDefaultScale * numLevelsAfterBoot,
	)

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: logN,
		LogQ: logQ,
		LogP: logP,
		LogDefaultScale: logDefaultScale,
		Xs: ring.Ternary{H: secretWeight},
	})

	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.RepackImagAsReal,
		LogSlots: logNumSlots,
		LevelQ: numLevelsAfterBoot + len(logQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		Levels: make([]int, len(logQiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {SlotsToCoeffsParameters.Levels[i] = 1}

	SCOREParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		LogSlots: logNumSlots,
		LevelQ: numLevelsAfterBoot + len(logQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		Levels: make([]int, len(logQiSlotsToCoeffs)),
	}
	for i := range SCOREParameters.Levels {SCOREParameters.Levels[i] = 1}

	bootParams_SPRU := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
	}

	bootParams_RSPRU := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SCOREParameters,
	}

	// Key generation & alike

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenSPRUKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	encoder := ckks.NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, pk)
	decryptor := rlwe.NewDecryptor(params, sk)

	evk_SPRU, _, _ := bootParams_SPRU.GenEvaluationKeys(sk)
	eval_SPRU, _ := bootstrapping.NewSlotsToCoeffsEvaluator(bootParams_SPRU, evk_SPRU)
	evk_RSPRU, _, _ := bootParams_RSPRU.GenEvaluationKeys(sk)
	eval_RSPRU, _ := bootstrapping.NewSCOREEvaluator(bootParams_RSPRU, evk_RSPRU)

	galoisElem := make([]uint64, logN+1) // For rotations and conjugation
	for i := range logN {galoisElem[i] = params.GaloisElementForRotation(1 << i)}
	galoisElem[logN] = params.GaloisElementForComplexConjugation()
	galoisKeys := kgen.GenGaloisKeysNew(galoisElem, sk)
	defaultEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, galoisKeys...))

	// Bootstrapping keys cs for SPRU

	baseRing := params.RingQ().AtLevel(0)
	sk_INTT := *sk.Value.Q.CopyNew()
	baseRing.INTT(sk_INTT, sk_INTT)
	baseRing.IMForm(sk_INTT, sk_INTT)

	csEncryptions_SPRU := make([]*rlwe.Ciphertext, 4*n)
	for u := range (4*n) {
		s := make([]complex128, N/2)
		for a := range (2*n) {
			for b := range secretWeight {
				for k := range Bover4n {
					i := 2*k*n*secretWeight + 2*b*n + a
					j := b*B + u*Bover4n + k
					s[extendedBitReverse(i, logNumSlots)] = complex(float64(sk_INTT.Coeffs[0][j]), 0)
				}
			}
		}
		sPoly := ckks.NewPlaintext(params, params.MaxLevel())
		sPoly.Scale = bootScale
		encoder.Encode(s, sPoly)
		cs, _ := encryptor.EncryptNew(sPoly)
		csEncryptions_SPRU[u] = cs
	}

	// Bootstrapping keys cs for R-SPRU

	csEncryptions_RSPRU := make([]*rlwe.Ciphertext, 2*n)
	for u := range (2*n) {
		s := make([]complex128, N/2)
		for a := range n {
			for b := range secretWeight {
				for k := range Bover2n {
					i := k*n*secretWeight + b*n + a
					j := b*B + u*Bover2n + k
					s[extendedBitReverse(i, logNumSlots)] = complex(float64(sk_INTT.Coeffs[0][j]), 0)
				}
			}
		}
		sPoly := ckks.NewPlaintext(params, params.MaxLevel())
		sPoly.Scale = bootScale
		encoder.Encode(s, sPoly)
		cs, _ := encryptor.EncryptNew(sPoly)
		csEncryptions_RSPRU[u] = cs
	}

	// Miscellaneous precomputations

	baseQ := int64(params.Q()[0])

	logExpScale := float64(logQBase[0]) - float64(logBootScale) - math.Log2(4*math.Pi)
	expScale := complex(math.Exp2(logExpScale / float64(secretWeight)), 0) // Called delta in paper
	
	// Test

	vecBeforeBoot := make([]complex128, N/2)
	vecBoot := make([]complex128, N/2)
	polyBeforeBoot := ckks.NewPlaintext(params, 0)

	var totalEncoding, totalCSProd, totalTrace, totalProdTree,
		totalSTC, totalSCORE, totalBoot, totalPrec float64
	var header string

	fmt.Println("╔══════════════════════╗")
	fmt.Println("║      SPRU tests      ║")
	fmt.Println("╚══════════════════════╝")
	fmt.Println()

	totalEncoding = 0
	totalCSProd = 0
	totalTrace = 0
	totalProdTree = 0
	totalSTC = 0
	totalBoot = 0
	totalPrec = 0

	header = fmt.Sprintf(
		"Run | %13s | %13s | %13s | %13s | %13s | %13s | %13s",
		"Encoding (s)",
		"CSProd (s)",
		"Trace (s)",
		"ProdTree (s)",
		"STC (s)",
		"Total (s)",
		"Precision",
	)
	fmt.Println(header)
	fmt.Println(strings.Repeat("─", len(header)))

	for run := range numBootRuns+1 {
		if run > 0 {fmt.Printf("%3d | ", run)} // Ignore the warm-up run
		for i := range n {vecBeforeBoot[i] = complex(sampling.RandFloat64(-1, 1), 0)}
		for i := n; i < len(vecBeforeBoot); i++ {vecBeforeBoot[i] = vecBeforeBoot[i%n]}
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)
		ctBeforeBoot, _ := encryptor.EncryptNew(polyBeforeBoot)

		// Encodings related to ciphertext coefficients

		t0 := time.Now()
	
		c0_INTT := *ctBeforeBoot.Value[0].CopyNew()
		c1_INTT := *ctBeforeBoot.Value[1].CopyNew()
		baseRing.INTT(c0_INTT, c0_INTT)
		baseRing.INTT(c1_INTT, c1_INTT)

		eEncodings := make([]*rlwe.Plaintext, 4*n)
		for u := range (4*n) {
			e := make([]complex128, N/2)
			for a := range (2*n) {
				l := a * N / (2 * n)
				for b := range secretWeight {
					for k := range Bover4n {
						i := 2*k*n*secretWeight + 2*b*n + a
						j := b*B + u*Bover4n + k
						var entry int64
						if j == 0 {
							entry = int64(c0_INTT.Coeffs[0][l] + c1_INTT.Coeffs[0][l])
						} else if j <= l {
							entry = int64(c1_INTT.Coeffs[0][l-j])
						} else {
							entry = baseQ - int64(c1_INTT.Coeffs[0][l-j+N])
						}
						e[extendedBitReverse(i, logNumSlots)] = expScale * psi(entry, baseQ)
					}
				}
			}
			ePoly := ckks.NewPlaintext(params, params.MaxLevel())
			ePoly.Scale = bootScale
			encoder.Encode(e, ePoly)
			eEncodings[u] = ePoly
		}

		// Multiplication with the cs

		t1 := time.Now()

		ctBoot := rlwe.NewCiphertext(params, 1, params.MaxLevel())
		for u := range (4*n) {
			term, _ := defaultEval.MulNew(csEncryptions_SPRU[u], eEncodings[u])
			if u == 0 {ctBoot = term} else {ctBoot, _ = defaultEval.AddNew(ctBoot, term)}
		}
		defaultEval.Rescale(ctBoot, ctBoot)

		// Trace

		t2 := time.Now()

		for i := N/4; i >= 2 * n * secretWeight; i /= 2 {
			ctRot, _ := defaultEval.RotateNew(ctBoot, i)
			ctBoot, _ = defaultEval.AddNew(ctBoot, ctRot)
		}

		// Product tree

		t3 := time.Now()

		for i := n * secretWeight; i >= 2 * n; i /= 2 {
			ctRot, _ := defaultEval.RotateNew(ctBoot, i)
			defaultEval.MulRelin(ctBoot, ctRot, ctBoot)
			defaultEval.Rescale(ctBoot, ctBoot)
		}

		// Imaginary part & SlotsToCoeffs

		t4 := time.Now()

		ctConj, _ := defaultEval.ConjugateNew(ctBoot)
		ctBoot, _ = defaultEval.SubNew(ctBoot, ctConj)
		eval_SPRU.Mul(ctBoot, -1i, ctBoot) // Division by imaginary unit
		ctBoot, _ = eval_SPRU.SlotsToCoeffs(ctBoot, nil)
		ctBoot.Scale = defaultScale

		// Computing the timings

		t5 := time.Now()

		timeEncoding := t1.Sub(t0).Seconds()
		timeCSProd := t2.Sub(t1).Seconds()
		timeTrace := t3.Sub(t2).Seconds()
		timeProdTree := t4.Sub(t3).Seconds()
		timeSTC := t5.Sub(t4).Seconds()
		timeBoot := t5.Sub(t0).Seconds()

		if run > 0 {
			totalEncoding += timeEncoding
			totalCSProd += timeCSProd
			totalTrace += timeTrace
			totalProdTree += timeProdTree
			totalSTC += timeSTC
			totalBoot += timeBoot

			encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
			stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
			prec := stats.AVGLog2Prec.Real
			totalPrec += prec

			fmt.Printf("%13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n",
				timeEncoding,
				timeCSProd,
				timeTrace,
				timeProdTree,
				timeSTC,
				timeBoot,
				prec,
			)
		}

		runtime.GC()
		debug.FreeOSMemory()
	}

	fmt.Println(strings.Repeat("─", len(header)))
	fmt.Printf("avg | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n\n",
		totalEncoding/float64(numBootRuns),
		totalCSProd/float64(numBootRuns),
		totalTrace/float64(numBootRuns),
		totalProdTree/float64(numBootRuns),
		totalSTC/float64(numBootRuns),
		totalBoot/float64(numBootRuns),
		totalPrec/float64(numBootRuns),
	)

	fmt.Println("╔══════════════════════╗")
	fmt.Println("║     R-SPRU tests     ║")
	fmt.Println("╚══════════════════════╝")
	fmt.Println()

	totalEncoding = 0
	totalCSProd = 0
	totalTrace = 0
	totalProdTree = 0
	totalSCORE = 0
	totalBoot = 0
	totalPrec = 0

	header = fmt.Sprintf(
		"Run | %13s | %13s | %13s | %13s | %13s | %13s | %13s",
		"Encoding (s)",
		"CSProd (s)",
		"Trace (s)",
		"ProdTree (s)",
		"SCORE (s)",
		"Total (s)",
		"Precision",
	)
	fmt.Println(header)
	fmt.Println(strings.Repeat("─", len(header)))

	for run := range numBootRuns+1 {
		if run > 0 {fmt.Printf("%3d | ", run)} // Ignore the warm-up run
		for i := range n {vecBeforeBoot[i] = complex(sampling.RandFloat64(-1, 1), 0)}
		for i := n; i < len(vecBeforeBoot); i++ {vecBeforeBoot[i] = vecBeforeBoot[i%n]}
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)
		ctBeforeBoot, _ := encryptor.EncryptNew(polyBeforeBoot)

		// Encodings related to ciphertext coefficients

		t0 := time.Now()
		
		c0_INTT := *ctBeforeBoot.Value[0].CopyNew()
		c1_INTT := *ctBeforeBoot.Value[1].CopyNew()
		baseRing.INTT(c0_INTT, c0_INTT)
		baseRing.INTT(c1_INTT, c1_INTT)

		eEncodings := make([]*rlwe.Plaintext, 2*n)
		for u := range (2*n) {
			e := make([]complex128, N/2)
			for a := range n {
				l := a * N / (2 * n)
				for b := range secretWeight {
					for k := range Bover2n {
						i := k*n*secretWeight + b*n + a
						j := b*B + u*Bover2n + k
						var entry int64
						if j == 0 {
							entry = int64(c0_INTT.Coeffs[0][l] + c1_INTT.Coeffs[0][l])
						} else if j <= l {
							entry = int64(c1_INTT.Coeffs[0][l-j])
						} else {
							entry = baseQ - int64(c1_INTT.Coeffs[0][l-j+N])
						}
						e[extendedBitReverse(i, logNumSlots)] = expScale * psi(entry, baseQ)
					}
				}
			}
			ePoly := ckks.NewPlaintext(params, params.MaxLevel())
			ePoly.Scale = bootScale
			encoder.Encode(e, ePoly)
			eEncodings[u] = ePoly
		}

		// Multiplication with the cs

		t1 := time.Now()

		ctBoot := rlwe.NewCiphertext(params, 1, params.MaxLevel())
		for u := range (2*n) {
			term, _ := defaultEval.MulNew(csEncryptions_RSPRU[u], eEncodings[u])
			if u == 0 {ctBoot = term} else {ctBoot, _ = defaultEval.AddNew(ctBoot, term)}
		}
		defaultEval.Rescale(ctBoot, ctBoot)

		// Trace

		t2 := time.Now()

		for i := N/4; i >= n * secretWeight; i /= 2 {
			ctRot, _ := defaultEval.RotateNew(ctBoot, i)
			ctBoot, _ = defaultEval.AddNew(ctBoot, ctRot)
		}

		// Product tree

		t3 := time.Now()

		for i := n * secretWeight / 2; i >= n; i /= 2 {
			ctRot, _ := defaultEval.RotateNew(ctBoot, i)
			defaultEval.MulRelin(ctBoot, ctRot, ctBoot)
			defaultEval.Rescale(ctBoot, ctBoot)
		}

		// Imaginary part & SCORE

		t4 := time.Now()

		ctConj, _ := defaultEval.ConjugateNew(ctBoot)
		ctBoot, _ = defaultEval.SubNew(ctBoot, ctConj)
		eval_RSPRU.Mul(ctBoot, -1i, ctBoot) // Division by imaginary unit
		ctBoot, _ = eval_RSPRU.SCORE(ctBoot)
		ctBoot.Scale = defaultScale

		// Computing the timings

		t5 := time.Now()

		timeEncoding := t1.Sub(t0).Seconds()
		timeCSProd := t2.Sub(t1).Seconds()
		timeTrace := t3.Sub(t2).Seconds()
		timeProdTree := t4.Sub(t3).Seconds()
		timeSCORE := t5.Sub(t4).Seconds()
		timeBoot := t5.Sub(t0).Seconds()

		if run > 0 {
			totalEncoding += timeEncoding
			totalCSProd += timeCSProd
			totalTrace += timeTrace
			totalProdTree += timeProdTree
			totalSCORE += timeSCORE
			totalBoot += timeBoot

			encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
			stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
			prec := stats.AVGLog2Prec.Real
			totalPrec += prec

			fmt.Printf("%13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n",
				timeEncoding,
				timeCSProd,
				timeTrace,
				timeProdTree,
				timeSCORE,
				timeBoot,
				prec,
			)
		}

		runtime.GC()
		debug.FreeOSMemory()
	}

	fmt.Println(strings.Repeat("─", len(header)))
	fmt.Printf("avg | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n\n",
		totalEncoding/float64(numBootRuns),
		totalCSProd/float64(numBootRuns),
		totalTrace/float64(numBootRuns),
		totalProdTree/float64(numBootRuns),
		totalSCORE/float64(numBootRuns),
		totalBoot/float64(numBootRuns),
		totalPrec/float64(numBootRuns),
	)
}

var bitReverseCache = make(map[[2]int]int)
func bitReverse(x, logN int) int {
	key := [2]int{x, logN}
	if val, ok := bitReverseCache[key]; ok {return val}
	res := 0
	for i := range logN {if (x>>i)&1 == 1 {res |= 1 << (logN - 1 - i)}}
	bitReverseCache[key] = res
	return res
}

var extendedBitReverseCache = make(map[[2]int]int)
func extendedBitReverse(x, logN int) int {
	key := [2]int{x, logN}
	if val, ok := extendedBitReverseCache[key]; ok {return val}
	N := 1 << logN
	res := (x/N)*N + bitReverse(x%N, logN)
	extendedBitReverseCache[key] = res
	return res
}

// Computing exp(2*pi*i * x/q)
func psi(x int64, q int64) complex128 {
	return cmplx.Rect(1, 2*math.Pi*float64(x)/float64(q))
}