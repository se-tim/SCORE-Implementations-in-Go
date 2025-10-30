# Modifications to Lattigo

**Upstream:** https://github.com/tuneinsight/lattigo  
**Base commit:** 5fc229b046cbba8ccbc00dd528af50e9de2307d5 (Sep 3, 2025)  
**Copied into:** `lattigo-main/`

## Summary of Changes

All modifications are summarized here; the original files remain unaltered elsewhere.

**`circuits/ckks/bootstrapping/evaluator.go`**:
- Added `SCOREMatrix` field to `Evaluator`.
- Added new constructors `NewSlotsToCoeffsEvaluator()` and `NewSCOREEvaluator()` to create evaluators for SlotsToCoeffs and SCORE.
- Updated `Evaluator.initialize()` to build SCORE matrices when the parameters are set.
- Added new method `Evaluator.SCORE()` for SCORE transformation.

**`circuits/ckks/dft/dft.go`**:
- Extended `MatrixLiteral` with `SCORE *SCOREOptions`.
- Added `SCOREOptions` structure.
- Added new methods `Evaluator.SCORENew()` and `Evaluator.SCORE()` implementing the SCORE transformation.
- Added SCORE-specific scaling logic in `MatrixLiteral.GenMatrices()`.

**`core/rlwe/keygenerator.go`**:
- Added SPRU key generation helpers:
  - `KeyGenerator.GenSPRUSecretKeyNew()`
  - `KeyGenerator.GenSPRUPairNew()`