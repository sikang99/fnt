package fnt

import (
	"fmt"
	"math"
	"testing"
)

func TestFHTLength(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	samples := randBlockF64(blockSize)

	fht := NewFHT(blockSize)

	defer func() {
		if r := recover(); r == nil {
			t.Fatal("failed to fail: " + r.(string))
		}
	}()
	fht.Execute(samples[1:], DIT, false)
	fht.Execute(samples[1:], DIF, false)
}

func TestFHTIdentity(t *testing.T) {
	bits := uint(15)
	blockSize := 1 << bits

	for _, d := range div {
		input := randBlockF64(blockSize)
		output := make([]float64, blockSize)
		copy(output, input)

		fht := NewFHT(blockSize)
		fht.Execute(output, d, false)
		fht.Execute(output, d, true)

		for i := range output {
			expected := input[i]
			received := output[i]
			diff := math.Abs(expected - received)
			if diff >= Tolerance {
				t.Fatalf("%f %f %f\n", diff, expected, received)
			}
		}
	}
}

func TestFHTConstant(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	fht := NewFHT(blockSize)

	for _, d := range div {
		input := make([]float64, blockSize)
		output := make([]float64, blockSize)
		for i := range input {
			input[i] = 1.0
			output[i] = 1.0
		}

		fht.Execute(output, d, false)

		if output[0] != float64(blockSize) {
			t.Fatalf("%f %f\n", float64(blockSize), output[0])
		}

		for i := 1; i < len(output)-1; i++ {
			if output[i] != 0.0 {
				t.Fatalf("%f %f\n", 0.0, output[i])
			}
		}
	}
}

func TestFHTDirect(t *testing.T) {
	bits := uint(7)
	blockSize := 1 << bits

	for _, d := range div {
		f0 := randBlockF64(blockSize)
		f1 := make([]float64, blockSize)
		copy(f1, f0)

		directHartleyTransform(f0)

		fht := NewFHT(blockSize)
		fht.Execute(f1, d, false)

		for i := range f0 {
			if math.Abs(f0[i]-f1[i]) > Tolerance {
				t.Fatalf("%f %f\n", f0[i], f1[i])
			}
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

func BenchmarkFHTRadix2DIT(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fht := NewFHT(blockSize)
	samples := randBlockF64(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fht.Execute(samples, DIT, false)
	}
}

func BenchmarkFHTRadix2DIF(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fht := NewFHT(blockSize)
	samples := randBlockF64(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fht.Execute(samples, DIF, false)
	}
}

func Example_hartley() {
	fht := NewFHT(8)

	samples := []float64{1, 1, 1, 1, 0, 0, 0, 0}
	fmt.Printf("%+0.3f\n", samples)

	// Transform to Hartley domain.
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
