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
	
	LogN := 14
	LogSlots := 13
	LogDefaultScale := 45

	q0 := []int{56}
	qiAfterBoot := []int{LogDefaultScale}
	qiSlotsToCoeffs := []int{39, 39, 39}
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}
	qiCoeffsToSlots := []int{56, 56, 56, 56} // Larger primes because of q-terms, more primes to make operation faster (is at higher modulus)
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
		LogSlots: LogSlots,
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
		LogSlots: LogSlots,
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

	fmt.Printf("logQP = %d\n", int(math.Round(float64(btpParams.BootstrappingParameters.LogQP()))))

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)
	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, _ := bootstrapping.NewEvaluator(btpParams, evk)

	// Encrypt

	myVector := make([]complex128, 1<<LogSlots)
	for i := range myVector {myVector[i] = sampling.RandComplex128(-1, 1)}
	myPoly := ckks.NewPlaintext(params, len(qiSlotsToCoeffs)) // Encrypt at lowest possible level
	encoder.Encode(myVector, myPoly)
	ct, _ := encryptor.EncryptNew(myPoly)
	encoder.Decode(decryptor.DecryptNew(ct), myVector) // Include encryption error in myVector

	// Bootstrap

	start := time.Now()

	ct, _ = eval.SlotsToCoeffs(ct, nil)
	ct, _, _ = eval.ScaleDown(ct) // Ensure the right gap between modulus and message
	ct, _ = eval.ModUp(ct)
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
	fmt.Printf("Bootstrapping time: %v\n", elapsed.Round(time.Millisecond))

	// Check

	myBootVector := make([]complex128, ct.Slots())
	encoder.Decode(decryptor.DecryptNew(ct), myBootVector)
	precStats := ckks.GetPrecisionStats(params, encoder, nil, myVector, myBootVector, 0, false)
	fmt.Println(precStats.String())

	fmt.Printf("eval.CoeffsToSlotsParameters.LogSlots = %d\n", eval.CoeffsToSlotsParameters.LogSlots)
}