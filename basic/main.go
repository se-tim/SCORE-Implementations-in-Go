package main

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/examples"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func main() {
	// Parameters
	params, _ := ckks.NewParametersFromLiteral(examples.CKKSComplexParamsN15QP881)

	// Keys & tools
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	encoder := ckks.NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, pk)
	decryptor := rlwe.NewDecryptor(params, sk)
	evaluator := ckks.NewEvaluator(params, evk)

	// Encode & encrypt
	slots := params.MaxSlots()
	values1 := make([]complex128, slots)
	values2 := make([]complex128, slots)
	values1[0] = 2
	values2[0] = 3

	pt1 := ckks.NewPlaintext(params, params.MaxLevel())
	pt2 := ckks.NewPlaintext(params, params.MaxLevel())
	encoder.Encode(values1, pt1)
	encoder.Encode(values2, pt2)

	ct1, _ := encryptor.EncryptNew(pt1)
	ct2, _ := encryptor.EncryptNew(pt2)

	// Homomorphic operations
	ctSum, _ := evaluator.AddNew(ct1, ct2)

	ctMul, _ := evaluator.MulNew(ct1, ct2)
	evaluator.Relinearize(ctMul, ctMul)
	evaluator.Rescale(ctMul, ctMul)

	// Decrypt & decode
	ptSum := decryptor.DecryptNew(ctSum)
	ptMul := decryptor.DecryptNew(ctMul)

	valuesSum := make([]complex128, slots)
	valuesMul := make([]complex128, slots)
	
	encoder.Decode(ptSum, valuesSum)
	encoder.Decode(ptMul, valuesMul)

	fmt.Printf("Addition: %.6f\n", real(valuesSum[0]))
	fmt.Printf("Addition level: %d\n", ctSum.Level())
	fmt.Printf("Multiplication: %.6f\n", real(valuesMul[0]))
	fmt.Printf("Multiplication level: %d\n", ctMul.Level())
	
	// Print some coefficients of the polynomial underlying ptSum
	for i := 0; i < 2; i++ {
		fmt.Printf("Coeff[%d] = %d\n", i, ptSum.Value.Coeffs[0][i])
	}

	// Dropping levels
	k := 0 // New level
	evaluator.DropLevel(ctMul, ctMul.Level()-k)
	ptMul = decryptor.DecryptNew(ctMul)
	encoder.Decode(ptMul, valuesMul)
	fmt.Printf("New multiplication: %.6f\n", real(valuesMul[0]))
	fmt.Printf("New multiplication level: %d\n", ctMul.Level())
}