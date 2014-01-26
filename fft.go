package fnt

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Stores transform size and pre-computed factors.
type FFT struct {
	n        int
	log      uint
	forward  []complex128
	backward []complex128
}

// Given transform size and direction, pre-computes all factors. Panics when n is not a
// power of 2.
func NewFFT(n int) (fft FFT) {
	// Check transform length is power of 2
	fft.n = n
	if fft.n&(fft.n-1) != 0 {
		panic(fmt.Sprintf("length must be power of 2: %d", fft.n))
	}

	// Calculate log base 2 of transform size
	fft.log = log2(uint(fft.n))

	// Allocate factor array
	n2 := fft.n >> 1
	fft.forward = make([]complex128, n2)
	fft.backward = make([]complex128, n2)

	// Pre-compute factors
	phi := math.Pi / float64(n2)
	for i := 0; i < n2; i++ {
		fft.forward[i] = cmplx.Rect(1, -phi*float64(i))
		fft.backward[i] = cmplx.Rect(1, phi*float64(i))
	}
	return
}

// Transforms vector f to Fourier space. The div parameter determines division
// kind. If norm is true, will divide the transformed vector f by the
// transform size. Panics when array length differs from expected length.
func (fft FFT) Execute(f []complex128, div DivisionKind, dir Direction, norm bool) {
	// Check given vector is expected length
	if len(f) != fft.n {
		panic("invalid transform length")
	}

	switch div {
	case DIT:
		fft.revBinPermute(f)
		if dir == Forward {
			fft.dit(f, fft.forward)
		} else {
			fft.dit(f, fft.backward)
		}
	case DIF:
		if dir == Forward {
			fft.dif(f, fft.forward)
		} else {
			fft.dif(f, fft.backward)
		}
		fft.revBinPermute(f)
	}

	if norm {
		fft.normalize(f)
	}
}

func (fft FFT) dit(f, factors []complex128) {
	for r := 0; r < fft.n; r += 2 {
		fft.sumDiff(&f[r], &f[r+1])
	}

	for ldm := uint(2); ldm <= fft.log; ldm++ {
		m := 1 << ldm
		mh := m >> 1
		stride := 1 << (fft.log - ldm)
		for r := 0; r < fft.n; r += m {
			for j, k := 0, 0; j < mh; j, k = j+1, k+stride {
				fft.sumDiffMultDIT(&f[r+j], &f[r+j+mh], factors[k])
			}
		}
	}
}

func (fft FFT) dif(f, factors []complex128) {
	for ldm := fft.log; ldm >= 2; ldm-- {
		m := 1 << ldm
		mh := m >> 1
		stride := 1 << (fft.log - ldm)
		for r := 0; r < fft.n; r += m {
			for j, k := 0, 0; j < mh; j, k = j+1, k+stride {
				fft.sumDiffMultDIF(&f[r+j], &f[r+j+mh], factors[k])
			}
		}

	}

	for r := 0; r < fft.n; r += 2 {
		fft.sumDiff(&f[r], &f[r+1])
	}
}

// Bit reversal permutation.
func (fft *FFT) revBinPermute(v []complex128) {
	var n, nh, r uint
	n = uint(len(v))
	nh = n >> 1

	for x := uint(1); x < nh; x++ {
		r += nh
		v[x], v[r] = v[r], v[x]
		x++

		for i := n; r&i == 0 && i > 0; {
			i >>= 1
			r ^= i
		}

		if r > x {
			v[x], v[r] = v[r], v[x]
			v[n-1-x], v[n-1-r] = v[n-1-r], v[n-1-x]
		}
	}
}

// Despite not modifying any state, passing a pointer avoids some overhead,
// which is substantial for these methods.
func (fft *FFT) sumDiff(a, b *complex128) {
	*a, *b = *a+*b, *a-*b
}

func (fft *FFT) sumDiffMultDIT(a, b *complex128, c complex128) {
	*a, *b = *a+*b*c, *a-*b*c
}

func (fft *FFT) sumDiffMultDIF(a, b *complex128, c complex128) {
	*a, *b = *a+*b, (*a-*b)*c
}

func (fft *FFT) normalize(f []complex128) {
	n := len(f)
	for i := range f {
		f[i] = f[i] / complex(float64(n), 0)
	}
}
