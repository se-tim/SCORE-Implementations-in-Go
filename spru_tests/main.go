package main

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
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
	
	LogN := 10
	LogSecretWeight := 4
	LogNumSlots := 3
	numLevelsAfterBoot := 1

	LogDefaultScale := 40
	LogBootScale := 58
	secretWeight := 1 << LogSecretWeight
	LogBaseQ := []int{55}
	LoqQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range LoqQiAfterBoot {LoqQiAfterBoot[i] = LogDefaultScale}
	LogQiSlotsToCoeffs := []int{39, 39, 39}
	LogQiProductTree := make([]int, LogSecretWeight)
	for i := range LogQiProductTree {LogQiProductTree[i] = LogBootScale}
	LogQiBootstrapKeyProduct := []int{LogBootScale}
	
	LogP := []int{61, 61, 61, 61, 61}
	LogQ := append(LogBaseQ, LoqQiAfterBoot...)
	LogQ = append(LogQ, LogQiSlotsToCoeffs...)
	LogQ = append(LogQ, LogQiProductTree...)
	LogQ = append(LogQ, LogQiBootstrapKeyProduct...)

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: LogN,
		LogQ: LogQ,
		LogP: LogP,
		LogDefaultScale: LogDefaultScale,
		Xs: ring.Ternary{H: secretWeight},
	})

	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicEncode,
		Format: dft.RepackImagAsReal,
		LogSlots: LogN - 1,
		LevelQ: params.MaxLevelQ(),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, 0),
	}
	for i := range CoeffsToSlotsParameters.Levels {
		CoeffsToSlotsParameters.Levels[i] = 1
	}

	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ: params.MaxLevelQ(),
		LogScale: LogBaseQ[0],
		Mod1Type: mod1.CosDiscrete,
		Mod1Degree: 0,
		DoubleAngle: 3,
		K: 16,
		LogMessageRatio: 24 - LogN,
		Mod1InvDegree: 0,
	}

	// Same parameters are used for SCORE
	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.RepackImagAsReal,
		LogSlots: LogNumSlots,
		LevelQ: numLevelsAfterBoot + len(LogQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: true,
		Levels: make([]int, len(LogQiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {
		SlotsToCoeffsParameters.Levels[i] = 1
	}

	btpParams := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		Mod1ParametersLiteral: Mod1ParametersLiteral,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		CircuitOrder: bootstrapping.ModUpThenEncode,
	}

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenSPRUKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	encoder := ckks.NewEncoder(params)
	// decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, err := bootstrapping.NewScoreEvaluator(btpParams, evk)
	if err != nil {panic(err)}

	rotElem := make([]uint64, LogN)
	for i := range LogN {rotElem[i] = params.GaloisElementForRotation(1 << i)}
	rotKeys := kgen.GenGaloisKeysNew(rotElem, sk)
	rotEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(rlk, rotKeys...))

	// Bootstrapping keys cs

	baseRing := params.RingQ().AtLevel(0)
	baseQ := uint64(params.Q()[0])

	sk_INTT := *sk.Value.Q.CopyNew()
	baseRing.INTT(sk_INTT, sk_INTT)
	baseRing.IMForm(sk_INTT, sk_INTT)

	n := 1 << LogNumSlots
	B := params.N() / secretWeight
	sVectors := make([][]complex128, 2*n)
	csEncryptions := make([]*rlwe.Ciphertext, 2*n)
	B_over_2n := B / (2 * n)

	for u := range (2 * n) {
		s := make([]complex128, 1<<(LogN-1))
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
		poly.Scale = rlwe.NewScale(math.Exp2(float64(LogBootScale)))
		encoder.Encode(s, poly)
		cs, _ := encryptor.EncryptNew(poly)
		sVectors[u] = s
		csEncryptions[u] = cs
	}

	// Define delta

	delta := math.Pow(
		(float64(uint64(1)<<LogBaseQ[0])/(4*math.Pi*math.Pow(2, float64(LogBootScale)))),
		1.0/float64(secretWeight),
	)

	// Encode & encrypt

	vec := make([]complex128, 1<<(LogN-1))
	for j := range vec {vec[j] = complex(float64(1+j%(1<<LogNumSlots)), 0)}
	poly := ckks.NewPlaintext(params, params.MaxLevel())
	encoder.Encode(vec, poly)
	ct, _ := encryptor.EncryptNew(poly)

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
		poly.Scale = rlwe.NewScale(math.Exp2(float64(LogBootScale)))
		encoder.Encode(e, poly)
		eEncodings[u] = poly
	}

	// Initial sum

	sum := rlwe.NewCiphertext(params, 1, params.MaxLevel())
	for u := range 2 * n {
		term, _ := rotEval.MulNew(csEncryptions[u], eEncodings[u])
		rotEval.Rescale(term, term)
		if u == 0 {sum = term} else {sum, _ = rotEval.AddNew(sum, term)}
	}

	// Trace

	for i := 1 << (LogN-1); i >= n * secretWeight; i /= 2 {
		ctRot, _ := rotEval.RotateNew(ct, i)
		ct, _ = rotEval.AddNew(ct, ctRot)
	}

	// Product tree

	for i := n * secretWeight / 2; i >= n; i /= 2 {
		ctRot, _ := rotEval.RotateNew(ct, i)
		rotEval.MulRelin(ct, ctRot, ct)
		rotEval.Rescale(ct, ct)
	}

	// SCORE

	ctConj, _ := rotEval.ConjugateNew(ct)
	ct, _ = rotEval.SubNew(ct, ctConj)
	ct, _ = eval.SCORE(ct)

	// Check

	decryptor := rlwe.NewDecryptor(params, sk)
	pt := decryptor.DecryptNew(ct)
	vecDec := make([]complex128, 1<<(LogN-1))
	encoder.Decode(pt, vecDec)
	fmt.Printf("Decrypted ct (first 10 values): %v\n", vecDec[:10])

	// Manual decryption check
	//
	// poly_INTT := *poly.Value.CopyNew()
	// baseRing.INTT(poly_INTT, poly_INTT)
	// baseRing.Reduce(poly_INTT, poly_INTT)

	// j := 0

	// fmt.Printf("poly_INTT.Coeffs[0][%d] = %d\n", j, poly_INTT.Coeffs[0][j]%baseQ)
	// 
	// res := (uint64(c0_INTT.Coeffs[0][j]) % baseQ)
	// for i := 0; i <= j; i++ {
	// 	term := mulMod(sk_INTT.Coeffs[0][i], c1_INTT.Coeffs[0][j-i], baseQ)
	// 	if term <= baseQ - res {res += term} else {res -= baseQ - term}
	// }
	// for i := j+1; i < params.N(); i++ {
	// 	term := mulMod(sk_INTT.Coeffs[0][i], c1_INTT.Coeffs[0][params.N()+j-i], baseQ)
	// 	if res >= term {res -= term} else {res += baseQ - term}
	// }
	// fmt.Printf("Estimated poly_INTT.Coeffs[0][%d] = %d\n", j, res)








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