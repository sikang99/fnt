package fnt

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"
)

func TestIdentity(t *testing.T) {
	bits := uint(15)
	blockSize := 1 << bits

	input := randBlock(blockSize)
	output := make([]float64, blockSize)
	copy(output, input)

	fht := NewFHT(bits)
	fht.Execute(output, false)
	fht.Execute(output, true)

	for i := range output {
		expected := input[i]
		received := output[i]
		diff := math.Abs(expected - received)
		if diff >= 1e-12 {
			t.Fatalf("%f %f %f\n", diff, expected, received)
		}
	}
}

func TestConstant(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	input := make([]float64, blockSize)
	output := make([]float64, blockSize)
	for i := range input {
		input[i] = 1.0
		output[i] = 1.0
	}

	fht := NewFHT(bits)
	fht.Execute(output, false)

	if output[0] != float64(blockSize) {
		t.Fatalf("%f %f\n", float64(blockSize), output[0])
	}

	for i := 1; i < len(output)-1; i++ {
		if output[i] != 0.0 {
			t.Fatalf("%f %f\n", 0.0, output[i])
		}
	}
}

func BenchmarkFHTRadix2(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fht := NewFHT(bits)
	samples := randBlock(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fht.Execute(samples, false)
	}
}

func ExampleFHT() {
	bits := uint(4)
	blockSize := 1 << bits

	samples := make([]float64, blockSize)
	for i := range samples {
		samples[i] = rand.Float64()
	}

	// NewFHT expects the log base 2 of the transform size.
	fht := NewFHT(bits)

	fmt.Printf("%0.3f\n", samples)
	fht.Execute(samples, false)
	fmt.Printf("%0.3f\n", samples)
}

func stepResponse(n int) []float64 {
	samples := make([]float64, n)
	for i := 0; i < n>>1; i++ {
		samples[i] = 1.0
	}
	return samples
}

func randBlock(n int) []float64 {
	samples := make([]float64, n)
	for i := range samples {
		samples[i] = rand.Float64()
	}
	return samples
}

func init() {
	rand.Seed(time.Now().UnixNano())
}
