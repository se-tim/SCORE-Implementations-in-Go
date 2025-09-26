package main

import (
	"flag"
	"fmt"
	"strings"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

type runResult struct {
    timeScaleDown float64
    timeModUp float64
    timeCTS float64
    timeEvalModRe float64
    timeEvalModIm float64
    timeSTC float64
    timeSCORE float64
    timeBoot float64
    prec float64
}

func main() {
	numBootRunsFlag := flag.Int("numBootRuns", 8, "number of bootstrap runs")
	LogNFlag := flag.Int("LogN", 14, "Log2 of ring degree N")
	levelsAfterBootFlag := flag.Int("levelsAfterBoot", 1, "levels after bootstrapping")
	flag.Parse()

	numBootRuns := *numBootRunsFlag
	LogN := *LogNFlag
	numLevelsAfterBoot := *levelsAfterBootFlag
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

	LogPQ := 0
	for _, q := range LogQ {LogPQ += q}
	for _, q := range LogP {LogPQ += q}

	println()
	fmt.Printf("LogN = %d, LogPQ = %d, LogModulusBeforeBoot = %d, LogModulusAfterBoot = %d\n\n",
		LogN,
		LogPQ,
		q0[0],
		q0[0] + LogDefaultScale * numLevelsAfterBoot,
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

	// Test

	var totalScaleDown, totalModUp, totalCTS, totalEvalModRe,
		totalEvalModIm, totalSTC, totalSCORE, totalBoot, totalPrec float64
	var header string
	var results []runResult
	var wg sync.WaitGroup

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

    results = make([]runResult, numBootRuns)
	wg = sync.WaitGroup{}

    for i := range numBootRuns {
        wg.Add(1)
        go func(i int) {
            defer wg.Done()

			kgen := rlwe.NewKeyGenerator(params)
			sk, pk := kgen.GenKeyPairNew()
			evk, _, _ := btpParams.GenEvaluationKeys(sk)
			encoder := ckks.NewEncoder(params)
			decryptor := rlwe.NewDecryptor(params, sk)
			encryptor := rlwe.NewEncryptor(params, pk)
			eval, _ := bootstrapping.NewEvaluator(btpParams, evk)

            vecBefore := make([]complex128, 1<<(LogN-1))
            for j := range vecBefore {
                vecBefore[j] = complex(sampling.RandFloat64(-1, 1), 0)
            }
			vecAfter := make([]complex128, 1<<(LogN-1))
            poly := ckks.NewPlaintext(params, 0)
            encoder.Encode(vecBefore, poly)
            ctBefore, _ := encryptor.EncryptNew(poly)
            encoder.Encode(vecBefore, poly)

            t0 := time.Now()
            ctBoot, _, _ := eval.ScaleDown(ctBefore)
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

            encoder.Decode(decryptor.DecryptNew(ctBoot), vecAfter)
            stats := ckks.GetPrecisionStats(params, encoder, nil, vecBefore, vecAfter, 0, false)
            prec := stats.AVGLog2Prec.Real

            results[i] = runResult{
                timeScaleDown: timeScaleDown,
                timeModUp: timeModUp,
                timeCTS: timeCTS,
                timeEvalModRe: timeEvalModRe,
                timeEvalModIm: timeEvalModIm,
                timeSTC: timeSTC,
                timeBoot: timeBoot,
                prec: prec,
            }

			fmt.Printf("Run %d/%d done...\n", i+1, numBootRuns)
        }(i)
    }

    wg.Wait()

	println()
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
        r := results[i]
        totalScaleDown += r.timeScaleDown
        totalModUp += r.timeModUp
        totalCTS += r.timeCTS
        totalEvalModRe += r.timeEvalModRe
        totalEvalModIm += r.timeEvalModIm
        totalSTC += r.timeSTC
        totalBoot += r.timeBoot
        totalPrec += r.prec
        fmt.Printf("%3d | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %13.3f | %8.1f bits\n",
            i+1, r.timeScaleDown, r.timeModUp, r.timeCTS, r.timeEvalModRe, r.timeEvalModIm, r.timeSTC, r.timeBoot, r.prec,
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

    results = make([]runResult, numBootRuns)
    wg = sync.WaitGroup{}

    for i := range numBootRuns {
        wg.Add(1)
        go func(i int) {
            defer wg.Done()

            kgen := rlwe.NewKeyGenerator(params)
			sk, pk := kgen.GenKeyPairNew()
			evk, _, _ := btpParams.GenEvaluationKeys(sk)
			encoder := ckks.NewEncoder(params)
			decryptor := rlwe.NewDecryptor(params, sk)
			encryptor := rlwe.NewEncryptor(params, pk)
			eval, _ := bootstrapping.NewEvaluator(btpParams, evk)

            vecBefore := make([]complex128, 1<<(LogN-1))
            for j := range vecBefore {
                vecBefore[j] = complex(sampling.RandFloat64(-1, 1), 0)
            }
			vecAfter := make([]complex128, 1<<(LogN-1))
            poly := ckks.NewPlaintext(params, 0)
            encoder.Encode(vecBefore, poly)
            ctBefore, _ := encryptor.EncryptNew(poly)
            encoder.Encode(vecBefore, poly)

            t0 := time.Now()
            ctBoot, _, _ := eval.ScaleDown(ctBefore)
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

            encoder.Decode(decryptor.DecryptNew(ctBoot), vecAfter)
            stats := ckks.GetPrecisionStats(params, encoder, nil, vecBefore, vecAfter, 0, false)
            prec := stats.AVGLog2Prec.Real

            results[i] = runResult{
                timeScaleDown: timeScaleDown,
                timeModUp: timeModUp,
                timeCTS: timeCTS,
                timeEvalModRe: timeEvalModRe,
                timeSCORE: timeSCORE,
                timeBoot: timeBoot,
                prec: prec,
            }

			fmt.Printf("Run %d/%d done...\n", i+1, numBootRuns)
        }(i)
    }

    wg.Wait()

	println()
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
        r := results[i]
        totalScaleDown += r.timeScaleDown
        totalModUp += r.timeModUp
        totalCTS += r.timeCTS
        totalEvalModRe += r.timeEvalModRe
        totalSCORE += r.timeSCORE
        totalBoot += r.timeBoot
        totalPrec += r.prec
        fmt.Printf("%3d | %13.3f | %13.3f | %13.3f | %13.3f | %13s | %13.3f | %13.3f | %8.1f bits\n",
            i+1, r.timeScaleDown, r.timeModUp, r.timeCTS, r.timeEvalModRe, "", r.timeSCORE, r.timeBoot, r.prec,
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