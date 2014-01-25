package fnt

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

func TestFFTLength(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	samples := randBlockC128(blockSize)

	fft := NewFFT(blockSize)

	defer func() {
		if r := recover(); r == nil {
			t.Fatal("failed to fail: " + r.(string))
		}
	}()
	fft.Execute(samples[1:], DIT, Forward, false)
	fft.Execute(samples[1:], DIF, Forward, false)
}

func TestFFTIdentity(t *testing.T) {
	bits := uint(14)
	blockSize := 1 << bits

	for _, d := range div {
		input := randBlockC128(blockSize)
		output := make([]complex128, blockSize)
		copy(output, input)

		fft := NewFFT(blockSize)
		fft.Execute(output, d, Forward, false)
		fft.Execute(output, d, Backward, true)

		for i := range output {
			expected := input[i]
			received := output[i]
			diff := cmplx.Abs(expected - received)
			if diff >= Tolerance {
				t.Fatalf("%f %f %f\n", diff, expected, received)
			}
		}
	}
}

func TestFFTConstant(t *testing.T) {
	bits := uint(4)
	blockSize := 1 << bits

	fft := NewFFT(blockSize)

	for d := range div {
		input := make([]complex128, blockSize)
		output := make([]complex128, blockSize)
		for i := range input {
			input[i] = 1.0
			output[i] = 1.0
		}

		fft.Execute(output, div[d], Forward, false)

		if output[0] != complex(float64(blockSize), 0) {
			t.Fatalf("%f %f\n", float64(blockSize), output[0])
		}

		for i := 1; i < len(output)-1; i++ {
			if output[i] != 0.0 {
				t.Fatalf("%f %f\n", 0.0, output[i])
			}
		}
	}
}

func TestFFTDirect(t *testing.T) {
	bits := uint(7)
	blockSize := 1 << bits

	for _, d := range div {
		f0 := randBlockC128(blockSize)
		f1 := make([]complex128, blockSize)
		copy(f1, f0)

		directFourierTransform(f0, -1.0)

		fft := NewFFT(blockSize)
		fft.Execute(f1, d, Forward, false)

		for i := range f0 {
			if cmplx.Abs(f0[i]-f1[i]) > Tolerance {
				t.Fatalf("%f %f\n", f0[i], f1[i])
			}
		}
	}
}

func directFourierTransform(f []complex128, sign float64) {
	n := len(f)
	h := make([]complex128, n)
	phi := sign * 2.0 * math.Pi / float64(n)
	for w := 0; w < n; w++ {
		var t complex128
		for k := 0; k < n; k++ {
			t += f[k] * cmplx.Rect(1, phi*float64(k)*float64(w))
		}
		h[w] = t
	}
	copy(f, h)
}

func BenchmarkFFTRadix2DIT(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fft := NewFFT(blockSize)
	samples := randBlockC128(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fft.Execute(samples, DIT, Forward, false)
	}
}

func BenchmarkFFTRadix2DIF(b *testing.B) {
	bits := uint(14)
	blockSize := 1 << bits

	fft := NewFFT(blockSize)
	samples := randBlockC128(blockSize)

	b.SetBytes(int64(blockSize))
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		fft.Execute(samples, DIF, Forward, false)
	}
}

func Example_fourier() {
	fft := NewFFT(8)

	samples := []complex128{1, 1, 1, 1, 0, 0, 0, 0}
	fmt.Printf("%+0.3f\n", samples)

	// Transform to Fourier domain.
	fft.Execute(samples, DIT, Forward, false)
	fmt.Printf("%+0.3f\n", samples)

	// Transform back to time-domain and normalize.
	fft.Execute(samples, DIT, Backward, true)
	fmt.Printf("%+0.3f\n", samples)
	// Output:
	// [(+1.000+0.000i) (+1.000+0.000i) (+1.000+0.000i) (+1.000+0.000i) (+0.000+0.000i) (+0.000+0.000i) (+0.000+0.000i) (+0.000+0.000i)]
	// [(+4.000+0.000i) (+1.000-2.414i) (+0.000+0.000i) (+1.000-0.414i) (+0.000+0.000i) (+1.000+0.414i) (+0.000+0.000i) (+1.000+2.414i)]
	// [(+1.000+0.000i) (+1.000+0.000i) (+1.000-0.000i) (+1.000-0.000i) (+0.000+0.000i) (+0.000-0.000i) (+0.000+0.000i) (+0.000+0.000i)]
}
