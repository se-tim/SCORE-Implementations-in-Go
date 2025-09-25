// CoeffsToSlots-first approach

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

	LogN := 12
	LogNumSlots := LogN - 1
	numLevelsAfterBoot := 10
	longTermSecretWeight := 192
	ephemeralSecretWeight := 32

	LogDefaultScale := 40
	q0 := []int{55}
	qiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range qiAfterBoot {qiAfterBoot[i] = LogDefaultScale}
	qiSlotsToCoeffs := []int{39, 39, 39}
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}
	qiCoeffsToSlots := []int{56, 56, 56, 56}
	LogP := []int{61, 61, 61, 61, 61}

	LogQ := append(q0, qiAfterBoot...)
	LogQ = append(LogQ, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: LogN,
		LogQ: LogQ,
		LogP: LogP,
		LogDefaultScale: LogDefaultScale,
		Xs: ring.Ternary{H: longTermSecretWeight},
	})

	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicEncode,
		Format: dft.RepackImagAsReal,
		LogSlots: LogNumSlots,
		LevelQ: params.MaxLevelQ(),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(qiCoeffsToSlots)),
	}
	for i := range CoeffsToSlotsParameters.Levels {CoeffsToSlotsParameters.Levels[i] = 1}

	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ: numLevelsAfterBoot + len(qiSlotsToCoeffs) + len(qiEvalMod),
		LogScale: qiEvalMod[0],
		Mod1Type: mod1.CosDiscrete,
		Mod1Degree: 30,
		DoubleAngle: 3,
		K: 16,
		LogMessageRatio: 24 - 2 * LogN + LogNumSlots,
		Mod1InvDegree: 0,
	}

	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.RepackImagAsReal,
		LogSlots: LogNumSlots,
		LevelQ: numLevelsAfterBoot + len(qiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(qiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {SlotsToCoeffsParameters.Levels[i] = 1}

	btpParams := bootstrapping.Parameters{
		ResidualParameters: params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral: Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight: ephemeralSecretWeight,
		CircuitOrder: bootstrapping.ModUpThenEncode,
	}

	// Key generation
	
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)
	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, _ := bootstrapping.NewEvaluator(btpParams, evk)
	fmt.Printf("LogN = %d, LogNumSlots = %d, LogPQ = %d, numLevelsAfterBoot = %d\n", LogN, LogNumSlots, int(math.Round(float64(btpParams.BootstrappingParameters.LogQP()))), numLevelsAfterBoot)

	// Test

	myVector := make([]complex128, 1<<LogNumSlots)
	for i := range myVector {myVector[i] = complex(sampling.RandFloat64(-1, 1), 0)}
	myPoly := ckks.NewPlaintext(params, 0)
	encoder.Encode(myVector, myPoly)
	ct, _ := encryptor.EncryptNew(myPoly)
	
	start := time.Now()
	ct, _, _ = eval.ScaleDown(ct)
	ct, _ = eval.ModUp(ct)
	ctReal, _, _ := eval.CoeffsToSlots(ct)
	ctReal, _ = eval.EvalMod(ctReal)
	ct, _ = eval.SCORE(ctReal)
	eval.Mul(ct, 1<<uint(LogN-1-LogNumSlots), ct) // Multiply by N / (2 * NumSlots)
	elapsed := time.Since(start)

	// Check

	myBootVector := make([]complex128, ct.Slots())
	encoder.Decode(decryptor.DecryptNew(ct), myBootVector)
	precStats := ckks.GetPrecisionStats(params, encoder, nil, myVector, myBootVector, 0, false)
	fmt.Println(precStats.String())
	fmt.Printf("myVector = [%f, %f, ...]\n", real(myVector[0]), real(myVector[1]))
	fmt.Printf("myBootVector = [%f, %f, ...]\n", real(myBootVector[0]), real(myBootVector[1]))
	fmt.Printf("Bootstrapping took %s\n", elapsed)
}