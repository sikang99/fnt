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

// Given log base 2 of the transform size, pre-computes all factors.
func NewFHT(n int) (fht FHT) {
	fht.n = n
	if fht.n&(fht.n-1) != 0 {
		panic(fmt.Sprintf("length must be power of 2: %d", fht.n))
	}

	for ; n > 1; fht.log, n = fht.log+1, n>>1 {
	}

	n2 := fht.n >> 1

	fht.factors = make([]float64, n2)

	phi := math.Pi / float64(n2)
	for idx := 0; idx < n2; idx += 2 {
		fht.factors[idx], fht.factors[idx+1] = math.Sincos(phi * float64((idx>>1)+1))
	}
	return
}

// Transforms vector f to Hartley space. If normalize is true, will divide the
// transformed vector f by the transform size. DIF transform is currently
// unimplemented.
func (fht FHT) Execute(f []float64, division DivisionKind, normalize bool) {
	if len(f) != fht.n {
		panic("invalid transform length")
	}

	switch division {
	case DIT:
		revBinPermute(f)
		fht.hartleyDIT(f, normalize)
	case DIF:
		// Unimplemented
		revBinPermute(f)
	}
}

func (fht FHT) hartleyDIT(f []float64, normalize bool) {
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
			for j, k := 1, m2-1; j < k; j, k = j+1, k-1 {
				s, c := fht.factors[fIdx], fht.factors[fIdx+1]

				sumDiffMult(&f[r+j+m2], &f[r+k+m2], s, c)
				sumDiff(&f[r+j], &f[r+j+m2])
				sumDiff(&f[r+k], &f[r+k+m2])

				fIdx += stride
			}

		}
	}

	if normalize {
		for i := range f {
			f[i] = f[i] / float64(n)
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
