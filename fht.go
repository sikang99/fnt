package fnt

import (
	"fmt"
	"math"
)

// Constants for selecting algorithm type.
type DivisionKind int

const (
	DIT DivisionKind = iota // Division-in-Time
	DIF                     // Division-in-Frequency
)

// Stores transform size and pre-computed factors.
type FHT struct {
	n       int
	log     uint
	factors []float64
}

// Given transform size, pre-computes all factors. Panics when n is not a
// power of 2.
func NewFHT(n int) (fht FHT) {
	// Check transform length is power of 2
	fht.n = n
	if fht.n&(fht.n-1) != 0 {
		panic(fmt.Sprintf("length must be power of 2: %d", fht.n))
	}

	// Calculate log base 2 of transform size
	for ; n > 1; fht.log, n = fht.log+1, n>>1 {
	}

	// Allocate factor array
	n2 := fht.n >> 1
	fht.factors = make([]float64, n2)

	// Pre-compute factors
	phi := math.Pi / float64(n2)
	for idx, i := 0, 1.0; idx < n2; idx, i = idx+2, i+1 {
		fht.factors[idx], fht.factors[idx+1] = math.Sincos(phi * i)
	}
	return
}

// Transforms vector f to Hartley space. The div parameter determines division
// kind. If norm is true, will divide the transformed vector f by the
// transform size. Panics when array length differs from expected length.
func (fht FHT) Execute(f []float64, div DivisionKind, norm bool) {
	// Check given vector is expected length
	if len(f) != fht.n {
		panic("invalid transform length")
	}

	switch div {
	case DIT:
		revBinPermute(f)
		fht.dit(f)
	case DIF:
		fht.dif(f)
		revBinPermute(f)
	}

	if norm {
		normalize(f)
	}
}

func (fht FHT) dit(f []float64) {
	n := 1 << fht.log

	for ldm := uint(1); ldm <= fht.log; ldm++ {
		m := 1 << ldm
		m2 := m >> 1
		m4 := m2 >> 1

		stride := 1 << (fht.log - ldm + 1)

		for r := 0; r < n; r += m {
			sumDiff(&f[r], &f[r+m2])

			if m4 != 0 {
				sumDiff(&f[r+m4], &f[r+m2+m4])
			}

			fIdx := stride - 2
			for j := 1; j < m4; j++ {
				s, c := fht.factors[fIdx], fht.factors[fIdx+1]

				sumDiffMult(&f[r+j+m2], &f[r+m2-j+m2], s, c)
				sumDiff(&f[r+j], &f[r+j+m2])
				sumDiff(&f[r+m2-j], &f[r+m2-j+m2])
				fIdx += stride
			}
		}
	}
}

func (fht FHT) dif(f []float64) {
	n := 1 << fht.log

	for ldm := fht.log; ldm >= 1; ldm-- {
		m := 1 << ldm
		m2 := m >> 1
		m4 := m2 >> 1

		stride := 1 << (fht.log - ldm + 1)

		for r := 0; r < n; r += m {
			for j := 0; j < m2; j++ {
				sumDiff(&f[r+j], &f[r+j+m2])
			}

			for j, fIdx := 1, stride-2; j < m4; j, fIdx = j+1, fIdx+stride {
				sumDiffMult(&f[r+m2+j], &f[r+m-j], fht.factors[fIdx], fht.factors[fIdx+1])
			}
		}
	}
}

// Bit reversal permutation.
func revBinPermute(v []float64) {
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

func sumDiff(a, b *float64) {
	*a, *b = *a+*b, *a-*b
}

func sumDiffMult(a, b *float64, s, c float64) {
	*a, *b = *a*c+*b*s, *a*s-*b*c
}

func normalize(f []float64) {
	n := len(f)
	for i := range f {
		f[i] = f[i] / float64(n)
	}
}
