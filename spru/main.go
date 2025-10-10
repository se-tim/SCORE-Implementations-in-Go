package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

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

// Compute exp(2*pi*i * x/q)
func psi(x int64, q int64) complex128 {
	return cmplx.Rect(1, 2 * math.Pi * float64(x) / float64(q))
}

func main() {
	// ==============
	//   Parameters
	// ==============

	runtime.GOMAXPROCS(1)
	
	logN := 14
	logSecretWeight := 6
	logNumSlots := logN - logSecretWeight - 2
	numLevelsAfterBoot := 1

	N := 1 << logN
	n := 1 << logNumSlots
	secretWeight := 1 << logSecretWeight
	B := N / secretWeight
	Bover4n := B / (4 * n)
	//Bover2n := 2 * Bover4n

	logDefaultScale := 40
	logBootScale := 60
	defaultScale := rlwe.NewScale(math.Exp2(float64(logDefaultScale)))
	bootScale := rlwe.NewScale(math.Exp2(float64(logBootScale)))

	logQBase := []int{53}
	logQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range logQiAfterBoot {logQiAfterBoot[i] = logDefaultScale}
	var logQiSlotsToCoeffs []int
	if logNumSlots <= 6 {logQiSlotsToCoeffs = []int{39}} else {logQiSlotsToCoeffs = []int{39, 39}}
	logQiRemaining := make([]int, 1+logSecretWeight+1) // Preparing STC, product tree, pt-ct products
	for i := range logQiRemaining {logQiRemaining[i] = logBootScale}
	
	logQ := append(logQBase, logQiAfterBoot...)
	logQ = append(logQ, logQiSlotsToCoeffs...)
	logQ = append(logQ, logQiRemaining...)

	logP := []int{61, 61, 61}

	logPQ := 0
	for _, q := range logQ {logPQ += q}
	for _, p := range logP {logPQ += p}
	
	println()
	fmt.Printf("logN = %d, logNumSlots = %d, logPQ = %d, logModulusBeforeBoot = %d, logModulusAfterBoot = %d\n",
		logN,
		logNumSlots,
		logPQ,
		logQBase[0],
		logQBase[0]+logDefaultScale*numLevelsAfterBoot,
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
		LogSlots: logNumSlots,
		LevelQ: numLevelsAfterBoot + len(logQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		Levels: make([]int, len(logQiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {SlotsToCoeffsParameters.Levels[i] = 1}

	btpParams := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
	}

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenSPRUKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, _ := bootstrapping.NewSlotsToCoeffsEvaluator(btpParams, evk)

	galoisElem := make([]uint64, logN+1)
	for i := range logN {galoisElem[i] = params.GaloisElementForRotation(1 << i)}
	galoisElem[logN] = params.GaloisElementForComplexConjugation()
	galoisKeys := kgen.GenGaloisKeysNew(galoisElem, sk)
	defaultEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, galoisKeys...))

	baseRing := params.RingQ().AtLevel(0)
	baseQ := int64(params.Q()[0])

	// Bootstrapping keys cs

	sk_INTT := *sk.Value.Q.CopyNew()
	baseRing.INTT(sk_INTT, sk_INTT)
	baseRing.IMForm(sk_INTT, sk_INTT)

	csEncryptions := make([]*rlwe.Ciphertext, 4*n)
	for u := range (4 * n) {
		s := make([]complex128, 1<<(logN-1))
		for a := range (2 * n) {
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
		csEncryptions[u] = cs
	}

	// Miscellaneous precomputations

	logExpScale := float64(logQBase[0]) - float64((1<<logSecretWeight)+logBootScale) - math.Log2(4*math.Pi)
	expScale := complex(math.Exp2(logExpScale / float64(secretWeight)), 0) // delta in paper

	maskVec := make([]complex128, N/2)
	for i := range (N/2) {if (i/n)%2 == 0 {maskVec[i] = 1} else {maskVec[i] = 1i}}
	maskPoly := ckks.NewPlaintext(params, params.MaxLevel())
	maskPoly.Scale = bootScale
	encoder.Encode(maskVec, maskPoly) // Required before applying SlotsToCoeffs

	// Encode & encrypt

	vecBeforeBoot := make([]complex128, N/2)
	for i := range n {vecBeforeBoot[i] = complex(sampling.RandFloat64(-1, 1), 0)}
	for i := n; i < len(vecBeforeBoot); i++ {vecBeforeBoot[i] = vecBeforeBoot[i % n]}

	polyBeforeBoot := ckks.NewPlaintext(params, 0)
	encoder.Encode(vecBeforeBoot, polyBeforeBoot)
	ctBeforeBoot, _ := encryptor.EncryptNew(polyBeforeBoot)

	// ========
	//   SPRU
	// ========

	// Encodings related to ciphertext coefficients

	t0 := time.Now()
	
	c0_INTT := *ctBeforeBoot.Value[0].CopyNew()
	c1_INTT := *ctBeforeBoot.Value[1].CopyNew()
	baseRing.INTT(c0_INTT, c0_INTT)
	baseRing.INTT(c1_INTT, c1_INTT)

	eEncodings := make([]*rlwe.Plaintext, 4*n)
	for u := range (4 * n) {
		e := make([]complex128, N/2)
		for a := range (2 * n) {
			k := a * N / (2 * n)
			for b := range secretWeight {
				for l := range Bover4n {
					i := 2*l*n*secretWeight + 2*b*n + a
					j := b*B + u*Bover4n + l
					var entry int64
					if j == 0 {
						entry = int64(c0_INTT.Coeffs[0][k] + c1_INTT.Coeffs[0][k])
					} else if j <= k {
						entry = int64(c1_INTT.Coeffs[0][k-j])
					} else {
						entry = baseQ - int64(c1_INTT.Coeffs[0][k-j+N])
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

	// Initial plaintext-ciphertext multiplications

	t1 := time.Now()

	ctBoot := rlwe.NewCiphertext(params, 1, params.MaxLevel())
	for u := range (4 * n) {
		term, _ := defaultEval.MulNew(csEncryptions[u], eEncodings[u])
		defaultEval.Rescale(term, term)
		if u == 0 {ctBoot = term} else {ctBoot, _ = defaultEval.AddNew(ctBoot, term)}
	}

	// Trace

	t2 := time.Now()

	for i := N/2; i >= 2 * n * secretWeight; i /= 2 {
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
	eval.Mul(ctBoot, -1i, ctBoot) // Division by imaginary unit
	defaultEval.Mul(ctBoot, maskPoly, ctBoot)
	defaultEval.Rescale(ctBoot, ctBoot)
	ctRotN, _ := defaultEval.RotateNew(ctBoot, n)
	ctBoot, _ = defaultEval.AddNew(ctBoot, ctRotN)

	ctBoot, err := eval.SlotsToCoeffs(ctBoot, nil)
	if err != nil {panic(err)}
	ctBoot.Scale = defaultScale

	// Computing timings

	t5 := time.Now()

	timeEncodings := t1.Sub(t0).Seconds()
	timeInitialSum := t2.Sub(t1).Seconds()
	timeTrace := t3.Sub(t2).Seconds()
	timeProductTree := t4.Sub(t3).Seconds()
	timeSCORE := t5.Sub(t4).Seconds()
	timeTotal := t5.Sub(t0).Seconds()

	// Check

	vecBoot := make([]complex128, N/2)
	encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
	stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
	prec := stats.AVGLog2Prec.Real

	fmt.Printf("Estimated precision after bootstrapping: %.1f bits\n", prec)

	fmt.Printf("Timings (seconds):\n")
	fmt.Printf("  Encodings:      %.3f\n", timeEncodings)
	fmt.Printf("  Initial Sum:    %.3f\n", timeInitialSum)
	fmt.Printf("  Trace:          %.3f\n", timeTrace)
	fmt.Printf("  Product Tree:   %.3f\n", timeProductTree)
	fmt.Printf("  SCORE:          %.3f\n", timeSCORE)
	fmt.Printf("  Total:          %.3f\n", timeTotal)
}