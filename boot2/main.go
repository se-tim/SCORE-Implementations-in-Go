package main

import (
	"fmt"
	"math"
	"runtime"
	"strings"
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
	runtime.GOMAXPROCS(runtime.NumCPU())

	// Parameters

	LogN := 16
	numLevelsAfterBoot := 7
	longTermSecretWeight := 192
	ephemeralSecretWeight := 32
	numBootRuns := 5

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
		LogSlots: LogN - 1,
		LevelQ: params.MaxLevelQ(),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(qiCoeffsToSlots)),
	}
	for i := range CoeffsToSlotsParameters.Levels {
		CoeffsToSlotsParameters.Levels[i] = 1
	}

	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ: numLevelsAfterBoot + len(qiSlotsToCoeffs) + len(qiEvalMod),
		LogScale: qiEvalMod[0],
		Mod1Type: mod1.CosDiscrete,
		Mod1Degree: 30,
		DoubleAngle: 3,
		K: 16,
		LogMessageRatio: 24 - LogN,
		Mod1InvDegree: 0,
	}

	// Same parameters are used for SCORE
	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type: dft.HomomorphicDecode,
		Format: dft.RepackImagAsReal,
		LogSlots: LogN - 1,
		LevelQ: numLevelsAfterBoot + len(qiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(qiSlotsToCoeffs)),
	}
	for i := range SlotsToCoeffsParameters.Levels {
		SlotsToCoeffsParameters.Levels[i] = 1
	}

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

	fmt.Printf("LogN = %d, LogPQ = %d, numLevelsAfterBoot = %d\n\n",
		LogN,
		int(math.Round(float64(btpParams.BootstrappingParameters.LogQP()))),
		numLevelsAfterBoot,
	)

	// Test

	vecBeforeBoot := make([]complex128, 1<<(LogN-1))
	vecBoot := make([]complex128, 1<<(LogN-1))
	polyBeforeBoot := ckks.NewPlaintext(params, 0)

	var totalScaleDown, totalModUp, totalCTS, totalEvalModRe,
		totalEvalModIm, totalSTC, totalSCORE, totalBoot, totalPrec float64
	var header string	

	fmt.Println("╔══════════════════════════════════════════════╗")
	fmt.Println("║     Results for conventional bootstrapping   ║")
	fmt.Println("╚══════════════════════════════════════════════╝")
	fmt.Println()

	totalScaleDown = 0.0
	totalModUp = 0.0
	totalCTS = 0.0
	totalEvalModRe = 0.0
	totalEvalModIm = 0.0
	totalSTC = 0.0
	totalBoot = 0.0
	totalPrec = 0.0

	header = fmt.Sprintf(
		"Run | %13s | %13s | %13s | %13s | %13s | %13s | %13s | %13s",
		"ScaleDown (s)",
		"ModUp (s)",
		"CTS (s)",
		"EvalModRe (s)",
		"EvalModIm (s)",
		"STC (s)",
		"Total (s)",
		"Precision",
	)
	fmt.Println(header)
	fmt.Println(strings.Repeat("─", len(header)))

	for i := range numBootRuns {
		fmt.Printf("%3d | ", i+1)
		for j := range vecBeforeBoot {
			vecBeforeBoot[j] = complex(sampling.RandFloat64(-1, 1), 0)
		}
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)
		ctBeforeBoot, _ := encryptor.EncryptNew(polyBeforeBoot)
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)

		// The following is equivalent to ctBoot = eval.Bootstrap(ctBeforeBoot)
		// but the individual steps are timed separately.
		t0 := time.Now()
		ctBoot, _, _ := eval.ScaleDown(ctBeforeBoot)
		t1 := time.Now()
		ctBoot, _ = eval.ModUp(ctBoot)
		t2 := time.Now()
		ctReal, ctImag, _ := eval.CoeffsToSlots(ctBoot)
		t3 := time.Now()
		ctReal, _ = eval.EvalMod(ctReal)
		t4 := time.Now()
		ctImag, _ = eval.EvalMod(ctImag)
		t5 := time.Now()
		ctBoot, _ = eval.SlotsToCoeffs(ctReal, ctImag)
		t6 := time.Now()

		timeScaleDown := float64(t1.Sub(t0)) / float64(time.Second)
		timeModUp := float64(t2.Sub(t1)) / float64(time.Second)
		timeCTS := float64(t3.Sub(t2)) / float64(time.Second)
		timeEvalModRe := float64(t4.Sub(t3)) / float64(time.Second)
		timeEvalModIm := float64(t5.Sub(t4)) / float64(time.Second)
		timeSTC := float64(t6.Sub(t5)) / float64(time.Second)
		timeBoot := float64(t6.Sub(t0)) / float64(time.Second)

		totalScaleDown += timeScaleDown
		totalModUp += timeModUp
		totalCTS += timeCTS
		totalEvalModRe += timeEvalModRe
		totalEvalModIm += timeEvalModIm
		totalSTC += timeSTC
		totalBoot += timeBoot

		encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
		stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
		prec := stats.AVGLog2Prec.Real
		totalPrec += prec

		fmt.Printf("%13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n",
			timeScaleDown, timeModUp, timeCTS, timeEvalModRe, timeEvalModIm, timeSTC, timeBoot, prec,
		)
	}

	fmt.Println(strings.Repeat("─", len(header)))
	fmt.Printf("avg | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n\n",
		totalScaleDown/float64(numBootRuns),
		totalModUp/float64(numBootRuns),
		totalCTS/float64(numBootRuns),
		totalEvalModRe/float64(numBootRuns),
		totalEvalModIm/float64(numBootRuns),
		totalSTC/float64(numBootRuns),
		totalBoot/float64(numBootRuns),
		totalPrec/float64(numBootRuns),
	)

	fmt.Println("╔══════════════════════════════════════════════╗")
	fmt.Println("║         Results for SCORE bootstrapping      ║")
	fmt.Println("╚══════════════════════════════════════════════╝")
	fmt.Println()

	totalScaleDown = 0.0
	totalModUp = 0.0
	totalCTS = 0.0
	totalEvalModRe = 0.0
	totalSCORE = 0.0
	totalBoot = 0.0
	totalPrec = 0.0

	header = fmt.Sprintf(
		"Run | %13s | %13s | %13s | %13s | %13s | %13s | %13s | %13s",
		"ScaleDown (s)",
		"ModUp (s)",
		"CTS (s)",
		"EvalModRe (s)",
		"",
		"SCORE (s)",
		"Total (s)",
		"Precision",
	)
	fmt.Println(header)
	fmt.Println(strings.Repeat("─", len(header)))

	for i := range numBootRuns {
		fmt.Printf("%3d | ", i+1)
		for j := range vecBeforeBoot {
			vecBeforeBoot[j] = complex(sampling.RandFloat64(-1, 1), 0)
		}
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)
		ctBeforeBoot, _ := encryptor.EncryptNew(polyBeforeBoot)
		encoder.Encode(vecBeforeBoot, polyBeforeBoot)

		t0 := time.Now()
		ctBoot, _, _ := eval.ScaleDown(ctBeforeBoot)
		t1 := time.Now()
		ctBoot, _ = eval.ModUp(ctBoot)
		t2 := time.Now()
		ctBoot, _, _ = eval.CoeffsToSlots(ctBoot)
		t3 := time.Now()
		ctBoot, _ = eval.EvalMod(ctBoot)
		t4 := time.Now()
		ctBoot, _ = eval.SCORE(ctBoot)
		t5 := time.Now()

		timeScaleDown := float64(t1.Sub(t0)) / float64(time.Second)
		timeModUp := float64(t2.Sub(t1)) / float64(time.Second)
		timeCTS := float64(t3.Sub(t2)) / float64(time.Second)
		timeEvalModRe := float64(t4.Sub(t3)) / float64(time.Second)
		timeSCORE := float64(t5.Sub(t4)) / float64(time.Second)
		timeBoot := float64(t5.Sub(t0)) / float64(time.Second)

		totalScaleDown += timeScaleDown
		totalModUp += timeModUp
		totalCTS += timeCTS
		totalEvalModRe += timeEvalModRe
		totalSCORE += timeSCORE
		totalBoot += timeBoot

		encoder.Decode(decryptor.DecryptNew(ctBoot), vecBoot)
		stats := ckks.GetPrecisionStats(params, encoder, nil, vecBeforeBoot, vecBoot, 0, false)
		prec := stats.AVGLog2Prec.Real
		totalPrec += prec

		fmt.Printf("%13.3f | %13.3f | %13.3f | %13.3f | %13s | %13.3f | %13.3f | %8.1f bits\n",
			timeScaleDown, timeModUp, timeCTS, timeEvalModRe, "", timeSCORE, timeBoot, prec,
		)
	}

	fmt.Println(strings.Repeat("─", len(header)))
	fmt.Printf("avg | %13.3f | %13.3f | %13.3f | %13.3f | %13s | %13.3f | %13.3f | %8.1f bits\n",
		totalScaleDown/float64(numBootRuns),
		totalModUp/float64(numBootRuns),
		totalCTS/float64(numBootRuns),
		totalEvalModRe/float64(numBootRuns),
		"",
		totalSCORE/float64(numBootRuns),
		totalBoot/float64(numBootRuns),
		totalPrec/float64(numBootRuns),
	)
}