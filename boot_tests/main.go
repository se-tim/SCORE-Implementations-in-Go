package main

import (
	"fmt"
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
	runtime.GOMAXPROCS(1)

	// Parameters

	LogN := 16
	numLevelsAfterBoot := 1
	longTermSecretWeight := 192
	ephemeralSecretWeight := 32
	numBootRuns := 3

	LogDefaultScale := 40
	LogBaseQ := []int{55}
	LogQiAfterBoot := make([]int, numLevelsAfterBoot)
	for i := range LogQiAfterBoot {LogQiAfterBoot[i] = LogDefaultScale}
	LogQiSlotsToCoeffs := []int{39, 39, 39}
	LogQiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}
	LogQiCoeffsToSlots := []int{56, 56, 56, 56}
	LogP := []int{61, 61, 61, 61, 61}

	LogQ := append(LogBaseQ, LogQiAfterBoot...)
	LogQ = append(LogQ, LogQiSlotsToCoeffs...)
	LogQ = append(LogQ, LogQiEvalMod...)
	LogQ = append(LogQ, LogQiCoeffsToSlots...)

	LogPQ := 0
	for _, q := range LogQ {LogPQ += q}
	for _, p := range LogP {LogPQ += p}

	println()
	fmt.Printf("LogN = %d, LogPQ = %d, LogModulusBeforeBoot = %d, LogModulusAfterBoot = %d\n\n",
		LogN,
		LogPQ,
		LogBaseQ[0],
		LogBaseQ[0] + LogDefaultScale * numLevelsAfterBoot,
	)

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
		Levels: make([]int, len(LogQiCoeffsToSlots)),
	}
	for i := range CoeffsToSlotsParameters.Levels {
		CoeffsToSlotsParameters.Levels[i] = 1
	}

	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ: numLevelsAfterBoot + len(LogQiSlotsToCoeffs) + len(LogQiEvalMod),
		LogScale: LogQiEvalMod[0],
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
		LevelQ: numLevelsAfterBoot + len(LogQiSlotsToCoeffs),
		LevelP: params.MaxLevelP(),
		LogBSGSRatio: 1,
		BitReversed: false,
		Levels: make([]int, len(LogQiSlotsToCoeffs)),
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
	encryptor := rlwe.NewEncryptor(params, pk)
	decryptor := rlwe.NewDecryptor(params, sk)
	evk, _, _ := btpParams.GenEvaluationKeys(sk)
	eval, _ := bootstrapping.NewEvaluator(btpParams, evk)

	// Test

	vecBeforeBoot := make([]complex128, 1<<(LogN-1))
	vecBoot := make([]complex128, 1<<(LogN-1))
	polyBeforeBoot := ckks.NewPlaintext(params, 0)

	var totalScaleDown, totalModUp, totalCTS, totalEvalModRe,
		totalEvalModIm, totalSTC, totalSCORE, totalBoot, totalPrec float64
	var header string	

	fmt.Println("╔══════════════════════╗")
	fmt.Println("║      BOOT tests      ║")
	fmt.Println("╚══════════════════════╝")
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

		timeScaleDown := t1.Sub(t0).Seconds()
		timeModUp := t2.Sub(t1).Seconds()
		timeCTS := t3.Sub(t2).Seconds()
		timeEvalModRe := t4.Sub(t3).Seconds()
		timeEvalModIm := t5.Sub(t4).Seconds()
		timeSTC := t6.Sub(t5).Seconds()
		timeBoot := t6.Sub(t0).Seconds()

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

	fmt.Println("╔══════════════════════╗")
	fmt.Println("║     R-BOOT tests     ║")
	fmt.Println("╚══════════════════════╝")
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

		timeScaleDown := t1.Sub(t0).Seconds()
		timeModUp := t2.Sub(t1).Seconds()
		timeCTS := t3.Sub(t2).Seconds()
		timeEvalModRe := t4.Sub(t3).Seconds()
		timeSCORE := t5.Sub(t4).Seconds()
		timeBoot := t5.Sub(t0).Seconds()

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