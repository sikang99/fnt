package fnt

import (
	"math/rand"
	"time"
)

// Constants for selecting algorithm type.
type DivisionKind int

const (
	DIT DivisionKind = iota // Division-in-Time
	DIF                     // Division-in-Frequency
)

var div = []DivisionKind{DIT, DIF}

type Direction int

const (
	Forward Direction = iota
	Backward
)

const Tolerance = 1e-12

func log2(n uint) (log uint) {
	for ; n > 1; log, n = log+1, n>>1 {
	}
	return
}

func randBlockF64(n int) []float64 {
	samples := make([]float64, n)
	for i := range samples {
		samples[i] = rand.Float64()
	}
	return samples
}

func randBlockC128(n int) []complex128 {
	samples := make([]complex128, n)
	for i := range samples {
		samples[i] = complex(rand.Float64(), rand.Float64())
	}
	return samples
}

func stepResponseF64(n int) []float64 {
	samples := make([]float64, n)
	for i := 0; i < n>>1; i++ {
		samples[i] = 1.0
	}
	return samples
}

func stepResponseC128(n int) []complex128 {
	samples := make([]complex128, n)
	for i := 0; i < n>>1; i++ {
		samples[i] = 1.0
	}
	return samples
}

func init() {
	rand.Seed(time.Now().UnixNano())
}
