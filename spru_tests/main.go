package main

import (
	"fmt"
	"math"
	"math/cmplx"

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
	
	logN := 13
	logSecretWeight := 6
	logNumSlots := 6
	numLevelsAfterBoot := 1

	N := 1 << logN
	n := 1 << logNumSlots
	secretWeight := 1 << logSecretWeight
	B := N / secretWeight
	B_over_2n := B / (2 * n)

	logDefaultScale := 40
	logBootScale := 60

	logBaseQ := []int{53}
	logQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range logQiAfterBoot {logQiAfterBoot[i] = logDefaultScale}
	var logQiSlotsToCoeffs []int
	if logNumSlots <= 2 {logQiSlotsToCoeffs = []int{39}} else {logQiSlotsToCoeffs = []int{39, 39, 39}}
	logQiProductTree := make([]int, logSecretWeight)
	for i := range logQiProductTree {logQiProductTree[i] = logBootScale}
	logQiBootstrapKeyProduct := []int{logBootScale}
	
	logQ := append(logBaseQ, logQiAfterBoot...)
	logQ = append(logQ, logQiSlotsToCoeffs...)
	logQ = append(logQ, logQiProductTree...)
	logQ = append(logQ, logQiBootstrapKeyProduct...)

	logP := []int{61, 61, 61, 61, 61}

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
	eval, _ := bootstrapping.NewSCOREEvaluator(btpParams, evk)

	galoisElem := make([]uint64, logN + 1)
	for i := range logN {galoisElem[i] = params.GaloisElementForRotation(1 << i)}
	galoisElem[logN] = params.GaloisElementForComplexConjugation()
	galoisKeys := kgen.GenGaloisKeysNew(galoisElem, sk)
	defaultEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, galoisKeys...))

	baseRing := params.RingQ().AtLevel(0)
	baseQ := int64(params.Q()[0])

	logExpScale := float64(logBaseQ[0]) - float64((1<<logSecretWeight)+logBootScale) - math.Log2(4*math.Pi)
	expScale := complex(math.Exp2(logExpScale / float64(secretWeight)), 0) // delta in paper

	// Bootstrapping keys cs

	sk_INTT := *sk.Value.Q.CopyNew()
	baseRing.INTT(sk_INTT, sk_INTT)
	baseRing.IMForm(sk_INTT, sk_INTT)

	csEncryptions := make([]*rlwe.Ciphertext, 2*n)
	for u := range (2 * n) {
		s := make([]complex128, 1<<(logN-1))
		for a := range n {
			for b := range secretWeight {
				for k := range (B / (2 * n)) {
					i := k*n*secretWeight + b*n + a
					j := b*B + u*B/(2*n) + k
					s[extendedBitReverse(i, logNumSlots)] = complex(float64(sk_INTT.Coeffs[0][j]), 0)
				}
			}
		}
		sPoly := ckks.NewPlaintext(params, params.MaxLevel())
		sPoly.Scale = rlwe.NewScale(math.Exp2(float64(logBootScale)))
		encoder.Encode(s, sPoly)
		cs, _ := encryptor.EncryptNew(sPoly)
		csEncryptions[u] = cs
	}

	// Encode & encrypt

	vecBeforeBoot := make([]complex128, N/2)
	for i := range n {vecBeforeBoot[i] = complex(sampling.RandFloat64(-1, 1), 0)}
	for i := n; i < len(vecBeforeBoot); i++ {vecBeforeBoot[i] = vecBeforeBoot[i % n]}

	polyBeforeBoot := ckks.NewPlaintext(params, params.MaxLevel())
	encoder.Encode(vecBeforeBoot, polyBeforeBoot)
	ct, _ := encryptor.EncryptNew(polyBeforeBoot)

	// ========
	//   SPRU
	// ========

	// Encodings related to ciphertext coefficients

	c0_INTT := *ct.Value[0].CopyNew()
	c1_INTT := *ct.Value[1].CopyNew()
	baseRing.INTT(c0_INTT, c0_INTT)
	baseRing.INTT(c1_INTT, c1_INTT)

	eEncodings := make([]*rlwe.Plaintext, 2*n)
	for u := range (2 * n) {
		e := make([]complex128, N/2)
		for a := range n {
			k := a * N / (2 * n)
			for b := range secretWeight {
				for l := range B_over_2n {
					i := l*n*secretWeight + b*n + a
					j := b*B + u*B_over_2n + l
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
		ePoly.Scale = rlwe.NewScale(math.Exp2(float64(logBootScale)))
		encoder.Encode(e, ePoly)
		eEncodings[u] = ePoly
	}

	// Initial sum

	ctBoot := rlwe.NewCiphertext(params, 1, params.MaxLevel())
	for u := range (2 * n) {
		term, _ := defaultEval.MulNew(csEncryptions[u], eEncodings[u])
		defaultEval.Rescale(term, term)
		if u == 0 {ctBoot = term} else {ctBoot, _ = defaultEval.AddNew(ctBoot, term)}
	}

	// Trace

	for i := N/2; i >= n * secretWeight; i /= 2 {
		ctRot, _ := defaultEval.RotateNew(ctBoot, i)
		ctBoot, _ = defaultEval.AddNew(ctBoot, ctRot)
	}

	// Product tree

	for i := n * secretWeight / 2; i >= n; i /= 2 {
		ctRot, _ := defaultEval.RotateNew(ctBoot, i)
		defaultEval.MulRelin(ctBoot, ctRot, ctBoot)
		defaultEval.Rescale(ctBoot, ctBoot)
	}

	// SCORE

	ctConj, _ := defaultEval.ConjugateNew(ctBoot)
	ctBoot, _ = defaultEval.SubNew(ctBoot, ctConj)
	eval.Mul(ctBoot, -1i, ctBoot) // Division by imaginary unit
	ctBoot, _ = eval.SCORE(ctBoot)
	ctBoot.Scale = rlwe.NewScale(math.Exp2(float64(logDefaultScale)))

	// Check

	vecBoot := make([]complex128, N/2)
	encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
	stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
	prec := stats.AVGLog2Prec.Real

	fmt.Printf("Estimated precision after bootstrapping: %.1f bits\n", prec)
}