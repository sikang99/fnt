package fnt

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"
)

func TestLength(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	samples := randBlock(blockSize)

	fht := NewFHT(blockSize)

	defer func() {
		if r := recover(); r == nil {
			t.Fatal("failed to fail: " + r.(string))
		}
	}()
	fht.Execute(samples[1:], DIT, false)
}

func TestIdentity(t *testing.T) {
	bits := uint(15)
	blockSize := 1 << bits

	input := randBlock(blockSize)
	output := make([]float64, blockSize)
	copy(output, input)

	fht := NewFHT(blockSize)
	fht.Execute(output, DIT, false)
	fht.Execute(output, DIT, true)

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

	fht := NewFHT(blockSize)
	fht.Execute(output, DIT, false)

	if output[0] != float64(blockSize) {
		t.Fatalf("%f %f\n", float64(blockSize), output[0])
	}

	for i := 1; i < len(output)-1; i++ {
		if output[i] != 0.0 {
			t.Fatalf("%f %f\n", 0.0, output[i])
		}
	}
}

func TestDirect(t *testing.T) {
	bits := uint(7)
	blockSize := 1 << bits

	f0 := randBlock(blockSize)
	f1 := make([]float64, blockSize)
	copy(f1, f0)

	directHartleyTransform(f0)

	fht := NewFHT(blockSize)
	fht.Execute(f1, DIT, false)

	for i := range f0 {
		if math.Abs(f0[i]-f1[i]) > 1e-12 {
			t.Fatalf("%f %f\n", f0[i], f1[i])
		}
	}
}

func directHartleyTransform(f []float64) {
	n := len(f)
	workspace := make([]float64, n)

	phi := 2.0 * math.Pi / float64(n)
	for w := 0; w < n; w++ {
		for k := range f {
			s, c := math.Sincos(phi * float64(k) * float64(w))
			workspace[w] += (c + s) * f[k]
		}
	}
	copy(f, workspace)
}

func BenchmarkFHTRadix2(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fht := NewFHT(blockSize)
	samples := randBlock(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fht.Execute(samples, DIT, false)
	}
}

func Example_hartley() {
	fht := NewFHT(8)

	samples := []float64{1, 1, 1, 1, 0, 0, 0, 0}
	fmt.Printf("%+0.3f\n", samples)

	// Transform to Hartley space.
	fht.Execute(samples, DIT, false)
	fmt.Printf("%+0.3f\n", samples)

	// Transform back to time-domain and normalize.
	fht.Execute(samples, DIT, true)
	fmt.Printf("%+0.3f\n", samples)
	// Output:
	// [+1.000 +1.000 +1.000 +1.000 +0.000 +0.000 +0.000 +0.000]
	// [+4.000 +3.414 +0.000 +1.414 +0.000 +0.586 +0.000 -1.414]
	// [+1.000 +1.000 +1.000 +1.000 +0.000 +0.000 +0.000 +0.000]
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
