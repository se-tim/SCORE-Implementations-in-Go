package main

import (
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func main() {
	// Parameters
	// 128-bit security with LogPQ <= 1200 for LogN = 16
	
	LogN := 16
	LogNumSlots := LogN - 1
	LogDefaultScale := 45
	numLevelsAfterBoot := 1

	q0 := []int{56}
	qiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range qiAfterBoot {
		qiAfterBoot[i] = LogDefaultScale
	}
	var qiSlotsToCoeffs, qiCoeffsToSlots []int
	switch {
	case LogNumSlots >= 14:
		qiSlotsToCoeffs = []int{39, 39, 39}
		qiCoeffsToSlots = []int{56, 56, 56}
	case LogNumSlots >= 9:
		qiSlotsToCoeffs = []int{39, 39}
		qiCoeffsToSlots = []int{56, 56, 56}
	case LogNumSlots >= 6:
		qiSlotsToCoeffs = []int{39}
		qiCoeffsToSlots = []int{56, 56}
	default:
		qiSlotsToCoeffs = []int{39}
		qiCoeffsToSlots = []int{56}
	}
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}
	LogP := []int{61, 61, 61, 61, 61}

	LogQ := append(q0, qiAfterBoot...)
	LogQ = append(LogQ, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	// Precomputations

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: LogN,
		LogQ: LogQ,
		LogP: LogP,
		LogDefaultScale: LogDefaultScale,
		Xs: ring.Ternary{H: 192},
	})

	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.Standard,
		LogSlots: LogNumSlots,
		LevelQ: len(qiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		Levels: make([]int, len(qiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {
		SlotsToCoeffsParameters.Levels[i] = 1
	}

	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicEncode,
		Format: dft.Standard,
		LogSlots: LogNumSlots,
		LevelQ: params.MaxLevelQ(),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		Levels: make([]int, len(qiCoeffsToSlots)),
	}
	for i := range CoeffsToSlotsParameters.Levels {
		CoeffsToSlotsParameters.Levels[i] = 1
	}

	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ: params.MaxLevel() - len(qiCoeffsToSlots),
		LogScale: qiEvalMod[0],
		Mod1Type: mod1.CosDiscrete,
		Mod1Degree: 30,
		DoubleAngle: 3,
		K: 16,
		LogMessageRatio: 10, // Gap between modulus and message
		Mod1InvDegree: 0,
	}

	btpParams := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral: Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight: 32,
		CircuitOrder: bootstrapping.DecodeThenModUp,
	}

	fmt.Printf("LogN = %d, LogNumSlots = %d, LogPQ = %d, numLevelsAfterBoot = %d\n", LogN, LogNumSlots, int(math.Round(float64(btpParams.BootstrappingParameters.LogQP()))), numLevelsAfterBoot)

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)
	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, _ := bootstrapping.NewEvaluator(btpParams, evk)

	// Encrypt

	myVector := make([]complex128, 1<<LogNumSlots)
	for i := range myVector {myVector[i] = sampling.RandComplex128(-1, 1)}
	myPoly := ckks.NewPlaintext(params, len(qiSlotsToCoeffs)) // Encrypt at lowest possible level
	encoder.Encode(myVector, myPoly)
	ct, _ := encryptor.EncryptNew(myPoly)
	encoder.Decode(decryptor.DecryptNew(ct), myVector) // Include encryption error in myVector

	// Bootstrap

	start := time.Now()

	ct, _ = eval.SlotsToCoeffs(ct, nil)
	eval.Evaluator.Mul(ct, 1<<uint(LogN-1-LogNumSlots), ct) // Multiply by N / (2 * NumSlots)
	ct, _, _ = eval.ScaleDown(ct) // Ensure the right gap between modulus and message
	ct, _ = eval.ModUp(ct) // Includes trace and division by N / ( 2 * NumSlots)
	ct, _, _ = eval.CoeffsToSlots(ct)
	ct_conj, _ := eval.Evaluator.ConjugateNew(ct)
	real, _ := eval.Evaluator.AddNew(ct, ct_conj)
	imag, _ := eval.Evaluator.SubNew(ct, ct_conj)
	eval.Evaluator.Mul(imag, -1i, imag)
	real, _ = eval.EvalModAndScale(real, 0.5)
	imag, _ = eval.EvalModAndScale(imag, 0.5)
	eval.Evaluator.Mul(imag, 1i, imag)
	eval.Evaluator.Add(real, imag, ct)

	elapsed := time.Since(start)
	fmt.Printf("Bootstrapping time: %v", elapsed.Round(time.Millisecond))

	// Check

	myBootVector := make([]complex128, ct.Slots())
	encoder.Decode(decryptor.DecryptNew(ct), myBootVector)
	precStats := ckks.GetPrecisionStats(params, encoder, nil, myVector, myBootVector, 0, false)
	fmt.Println(precStats.String())
}