package main

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func mulMod(a, b, m uint64) uint64 {
	var A, B, M, P big.Int
	A.SetUint64(a)
	B.SetUint64(b)
	M.SetUint64(m)
	P.Mul(&A, &B)
	P.Mod(&P, &M)
	return P.Uint64()
}

func myMod(x, q uint64) int64 {
	r := int64(x % q)
	halfQ := int64(q / 2)
	if r > halfQ {
		r -= int64(q)
	}
	return r
}

// bitReverse returns the bit-reversal of x with respect to log2(n) bits.
// n must be a power of two.
func bitReverse(x, n int) int {
	logN := 0
	for tmp := n; tmp > 1; tmp >>= 1 {
		logN++
	}
	res := 0
	for i := 0; i < logN; i++ {
		if (x>>i)&1 == 1 {
			res |= 1 << (logN - 1 - i)
		}
	}
	return res
}

func extendedBitReverse(x, n int) int {
	return (x/n)*n + bitReverse(x%n, n)
}

func bitReverseVector(vec []complex128, n int) []complex128 {
	res := make([]complex128, len(vec))
	for i := 0; i < len(vec); i++ {
		j := extendedBitReverse(i, n)
		res[j] = vec[i]
	}
	return res
}

// psi computes y * exp(2*pi*i * x/q) with high precision using math/big.
func psi(x int64, y float64, q int64) complex128 {
	prec := 256
	twoPi := new(big.Float).SetPrec(uint(prec)).SetFloat64(2 * math.Pi)
	X := new(big.Float).SetPrec(uint(prec)).SetInt64(x)
	Q := new(big.Float).SetPrec(uint(prec)).SetInt64(q)

	// theta = 2*pi*x/q
	theta := new(big.Float).Quo(new(big.Float).Mul(twoPi, X), Q)
	theta64, _ := theta.Float64()

	return cmplx.Rect(y, theta64)
}

func main() {
	// ==============
	//   Parameters
	// ==============
	
	logN := 12
	logSecretWeight := 4
	logNumSlots := 2
	numLevelsAfterBoot := 3

	logDefaultScale := 40
	logBootScale := 58
	secretWeight := 1 << logSecretWeight
	logBaseQ := []int{55}
	loqQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range loqQiAfterBoot {loqQiAfterBoot[i] = logDefaultScale}
	// logQiSlotsToCoeffs := []int{39, 39, 39}
	logQiSlotsToCoeffs := []int{39}
	logQiProductTree := make([]int, logSecretWeight)
	for i := range logQiProductTree {logQiProductTree[i] = logBootScale}
	logQiBootstrapKeyProduct := []int{logBootScale}
	
	logP := []int{61, 61, 61, 61, 61}
	logQ := append(logBaseQ, loqQiAfterBoot...)
	logQ = append(logQ, logQiSlotsToCoeffs...)
	logQ = append(logQ, logQiProductTree...)
	logQ = append(logQ, logQiBootstrapKeyProduct...)

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: logN,
		LogQ: logQ,
		LogP: logP,
		LogDefaultScale: logDefaultScale,
		Xs: ring.Ternary{H: secretWeight},
	})

	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.Standard,
		LogSlots: logNumSlots,
		LevelQ: numLevelsAfterBoot + len(logQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(logQiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {
		SlotsToCoeffsParameters.Levels[i] = 1
	}

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

	rotElem := make([]uint64, logN + 1)
	for i := range logN {rotElem[i] = params.GaloisElementForRotation(1 << i)}
	rotElem[logN] = params.GaloisElementForComplexConjugation()
	rotKeys := kgen.GenGaloisKeysNew(rotElem, sk)
	rotEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, rotKeys...))

	// Bootstrapping keys cs

	baseRing := params.RingQ().AtLevel(0)
	baseQ := uint64(params.Q()[0])

	sk_INTT := *sk.Value.Q.CopyNew()
	baseRing.INTT(sk_INTT, sk_INTT)
	baseRing.IMForm(sk_INTT, sk_INTT)

	n := 1 << logNumSlots
	B := params.N() / secretWeight
	sVectors := make([][]complex128, 2*n)
	csEncryptions := make([]*rlwe.Ciphertext, 2*n)
	B_over_2n := B / (2 * n)

	for u := range (2 * n) {
		s := make([]complex128, 1<<(logN-1))
		for a := range n {
			for b := range secretWeight {
				for k := range B_over_2n {
					i := k*n*secretWeight + b*n + a
					j := b*B + u*B/(2*n) + k
					s[i] = complex(float64(sk_INTT.Coeffs[0][j]), 0)
				}
			}
		}
		s = bitReverseVector(s, n)
		poly := ckks.NewPlaintext(params, params.MaxLevel())
		poly.Scale = rlwe.NewScale(math.Exp2(float64(logBootScale)))
		encoder.Encode(s, poly)
		cs, _ := encryptor.EncryptNew(poly)
		sVectors[u] = s
		csEncryptions[u] = cs
	}

	// Define delta

	delta := math.Pow(
		(math.Pow(2, float64(logBaseQ[0]))/(4*math.Pi*math.Pow(2, 16+float64(logBootScale)))),
		1.0/float64(secretWeight),
	)

	// Encode & encrypt

	vec := make([]complex128, 1<<(logN-1))
	for k := range vec {vec[k] = complex(float64(k%n), 0)}
	// Fill base pattern with random reals in [-1, 1]
	// for i := 0; i < 1 << logNumSlots; i++ {
	// 	vec[i] = complex(rand.Float64()*2-1, 0)
	// }
	// // Repeat the base pattern across the rest of the vector
	// for i := 1 << logNumSlots; i < 1 << (logN-1); i++ {
	// 	vec[i] = vec[i % (1 << logNumSlots)]
	// }
	poly := ckks.NewPlaintext(params, params.MaxLevel())
	encoder.Encode(vec, poly)
	ct, _ := encryptor.EncryptNew(poly)

	polyDec := decryptor.DecryptNew(ct)
	polyDec.Scale = rlwe.NewScale(math.Exp2(float64(logDefaultScale)))
	vecDecCheck := make([]complex128, params.N()/2)
	encoder.Decode(polyDec, vecDecCheck)
	fmt.Printf("First 4 values of decoded polyDec: %v\n", vecDecCheck[:6])

	// ========
	//   SPRU
	// ========

	// Encodings related to ciphertext coefficients (TODO: merge loops)

	c0_INTT := *ct.Value[0].CopyNew()
	c1_INTT := *ct.Value[1].CopyNew()
	baseRing.INTT(c0_INTT, c0_INTT)
	baseRing.INTT(c1_INTT, c1_INTT)

	cVectors := make([][]int64, n)
	N := params.N()
	for a := range n {
		c_a := make([]int64, N)
		j := a * N / (2 * n)
		c_a[0] = myMod(c0_INTT.Coeffs[0][j] + c1_INTT.Coeffs[0][j], baseQ)
		for i := 1; i <= j; i++ {
			c_a[i] = myMod(c1_INTT.Coeffs[0][j-i], baseQ)
		}
		for i := j+1; i < N; i++ {
			c_a[i] = myMod(baseQ-c1_INTT.Coeffs[0][j-i+N], baseQ)
		}
		cVectors[a] = c_a
	}
	eVectors := make([][]complex128, 2*n)
	eEncodings := make([]*rlwe.Plaintext, 2*n)
	for u := range 2 * n {
		e := make([]complex128, params.N()/2)
		for a := range n {
			for b := range secretWeight {
				for k := range B_over_2n {
					i := k*n*secretWeight + b*n + a
					j := b*B + u*B_over_2n + k
					e[i] = complex(delta, 0) * psi(cVectors[a][j], 1.0, int64(baseQ))
				}
			}
		}
		e = bitReverseVector(e, n)
		poly := ckks.NewPlaintext(params, params.MaxLevel())
		poly.Scale = rlwe.NewScale(math.Exp2(float64(logBootScale)))
		encoder.Encode(e, poly)
		eVectors[u] = e
		eEncodings[u] = poly
	}

	// Initial sum

	ctBoot := rlwe.NewCiphertext(params, 1, params.MaxLevel())
	for u := range (2 * n) {
		term, _ := rotEval.MulNew(csEncryptions[u], eEncodings[u])
		rotEval.Rescale(term, term)
		if u == 0 {ctBoot = term} else {ctBoot, _ = rotEval.AddNew(ctBoot, term)}
	}

	// Trace

	for i := 1 << (logN-1); i >= n * secretWeight; i /= 2 {
		ctRot, _ := rotEval.RotateNew(ctBoot, i)
		ctBoot, _ = rotEval.AddNew(ctBoot, ctRot)
	}

	// Product tree

	for i := n * secretWeight / 2; i >= n; i /= 2 {
		ctRot, _ := rotEval.RotateNew(ctBoot, i)
		rotEval.MulRelin(ctBoot, ctRot, ctBoot)
		rotEval.Rescale(ctBoot, ctBoot)
	}

	// SCORE

	ctConj, _ := rotEval.ConjugateNew(ctBoot)
	ctBoot, _ = rotEval.SubNew(ctBoot, ctConj)
	eval.Mul(ctBoot, -1i, ctBoot) // Division by imaginary unit

	// J := bitReverse(j,n) * params.N() / (2 * n)

	// for j := 0; j < 3+1<<logNumSlots; j++ {
	// 	fmt.Printf("Cipher entry after taking twice imaginary part [%d]: %v\n", j, vecDecCheck[j])
	// }

	ctBoot0, _ := eval.SCORE(ctBoot)
	// ctBootDouble, _ := rotEval.MulNew(ctBoot, 2)
	// ctBoot1, _ := eval.SCORE(ctBootDouble)

	polyDec = decryptor.DecryptNew(ctBoot0)
	polyDec.Scale = rlwe.NewScale(math.Exp2(float64(logDefaultScale)))
	vecDecCheck = make([]complex128, params.N()/2)
	encoder.Decode(polyDec, vecDecCheck)
	fmt.Printf("First 4 values of decoded polyDec: %v\n", vecDecCheck[:6])

	stats := ckks.GetPrecisionStats(params, encoder, nil, vec, vecDecCheck, 0, false)
	prec := stats.AVGLog2Prec.Real

	fmt.Printf("Estimated precision after bootstrapping: %8.1f bits\n", prec)

	// Check

	higherRing := params.RingQ().AtLevel(2)

	for j := 0; j < 2*n; j++ {
		J := j * params.N() / (2 * n)
		polyBoot0 := decryptor.DecryptNew(ctBoot0)
		polyBoot0_INTT := *polyBoot0.Value.CopyNew()
		higherRing.INTT(polyBoot0_INTT, polyBoot0_INTT)
		fmt.Printf("boot-plaintext polynomial coefficient [%d]: %d\n", J, polyBoot0_INTT.Coeffs[0][J])
	}

	// polyBoot1 := decryptor.DecryptNew(ctBoot1)
	// polyBoot1_INTT := *polyBoot1.Value.CopyNew()
	// higherRing.INTT(polyBoot1_INTT, polyBoot1_INTT)
	// fmt.Printf("Some boot-plaintext polynomial coefficient (double input): %d\n", polyBoot1_INTT.Coeffs[0][J])


	fmt.Printf("Level of ctBoot: %d\n", ctBoot.Level())

	// vecDec := make([]complex128, n)
	// encoder.Decode(polyBoot, vecDec)
	// fmt.Printf("First slot of decrypted ct: %v\n", vecDec[0])


	// ==================
	//   Plaintext SPRU
	// ==================

	// Initial sum
	plainSum := make([]complex128, params.N()/2)
	for u := range (2 * n) {
		for i := range plainSum {
			plainSum[i] += sVectors[u][i] * eVectors[u][i]
		}
	}

	// Trace
	for i := 1 << (logN-1); i >= n*secretWeight; i /= 2 {
		for j := range plainSum {
			plainSum[j] += plainSum[(j+i)%(params.N()/2)]
		}
	}

	// Product tree
	for i := n * secretWeight / 2; i >= n; i /= 2 {
		for j := range plainSum {
			plainSum[j] *= plainSum[(j+i)%(params.N()/2)]
		}
	}

	// Take imaginary part
	for i := range plainSum {
		plainSum[i] -= cmplx.Conj(plainSum[i])
		plainSum[i] *= complex(0, -1)
	}

	// for j := 0; j < 1+1<<logNumSlots; j++ {
	// 	fmt.Printf("Some plain entry after taking twice imaginary part [%d]: %v\n", j, plainSum[j])
	// }







	// Get some plaintext polynomial coefficient
	
	for j := 0; j < 2*n; j++ {
		J := j * params.N() / (2 * n)
		poly_INTT := *poly.Value.CopyNew()
		baseRing.INTT(poly_INTT, poly_INTT)
		baseRing.Reduce(poly_INTT, poly_INTT)
		fmt.Printf("plaintext polynomial coefficient [%d]: %d\n", J, poly_INTT.Coeffs[0][J]%baseQ)
	}

	// SPRU in plaintext







	// Manual decryption
	
	// res := (uint64(c0_INTT.Coeffs[0][J]) % baseQ)
	// for i := 0; i <= J; i++ {
	// 	term := mulMod(sk_INTT.Coeffs[0][i], c1_INTT.Coeffs[0][J-i], baseQ)
	// 	if term <= baseQ - res {res += term} else {res -= baseQ - term}
	// }
	// for i := J+1; i < params.N(); i++ {
	// 	term := mulMod(sk_INTT.Coeffs[0][i], c1_INTT.Coeffs[0][params.N()+J-i], baseQ)
	// 	if res >= term {res -= term} else {res += baseQ - term}
	// }
	// fmt.Printf("Estimated poly_INTT.Coeffs[0][%d] = %d\n", J, res)

	// Trace

	// for i := 1 << (LogNumSlots-1); i >= 1; i /= 2 {
	// 	ctRot, _ := rotEval.RotateNew(ct, i)
	// 	ct, _ = rotEval.AddNew(ct, ctRot)
	// }

	// pt := decryptor.DecryptNew(ct)
	// vecDec := make([]complex128, 1<<LogNumSlots)
	// encoder.Decode(pt, vecDec)
	// fmt.Printf("Decrypted ct: %v\n", vecDec)






	// Product

	// for i := 1 << (LogNumSlots-1); i >= 1<<(LogNumSlots-2); i /= 2 {
	// 	ctRot, _ := rotEval.RotateNew(ct, i)
	// 	rotEval.MulRelin(ct, ctRot, ct)
	// 	rotEval.Rescale(ct, ct)
	// }

	// pt := decryptor.DecryptNew(ct)
	// vecDec := make([]complex128, 1<<LogNumSlots)
	// encoder.Decode(pt, vecDec)
	// fmt.Printf("Decrypted ct: %v\n", vecDec)








	// ctSum, _ := eval.AddNew(ct, ctRot)
	// decryptor := rlwe.NewDecryptor(params, sk)
	// ptSum := decryptor.DecryptNew(ctSum)
	// vecSum := make([]complex128, 1<<(LogNumSlots-1))
	// encoder.Decode(ptSum, vecSum)
	// fmt.Printf("Decrypted ctSum (first 10 values): %v\n", vecSum[:10])

	// decryptor := rlwe.NewDecryptor(params, sk)

	// vecBeforeBoot := make([]complex128, 1<<(LogN-1))
	// for j := range vecBeforeBoot {
	// 	vecBeforeBoot[j] = complex(sampling.RandFloat64(-1, 1), 0)
	// }
	// fmt.Printf("vecBeforeBoot[0] = %v\n", vecBeforeBoot[0])
	// fmt.Printf("vecBeforeBoot[1] = %v\n", vecBeforeBoot[1])
	// encoder.Encode(vecBeforeBoot, poly)
	// ct_test, _ := encryptor.EncryptNew(poly)
	// encoder.Decode(decryptor.DecryptNew(ct_test), vecBeforeBoot)

	// fmt.Printf("vecBeforeBoot[0] = %v\n", vecBeforeBoot[0])
	// fmt.Printf("vecBeforeBoot[1] = %v\n", vecBeforeBoot[1])

	//decryptor := rlwe.NewDecryptor(params, sk)
	//pt := decryptor.DecryptNew(csList[5])

	// vecBeforeBoot := make([]complex128, 1<<(LogN-1))
	// for j := range vecBeforeBoot {
	// 		vecBeforeBoot[j] = complex(sampling.RandFloat64(-1, 1), 0)
	// 	}
	// poly = ckks.NewPlaintext(params, 3)
	// cs, _ = encryptor.EncryptNew(poly)
	// encoder.Decode(pt, vecBoot)
}