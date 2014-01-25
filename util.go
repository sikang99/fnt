package fnt

// Constants for selecting algorithm type.
type DivisionKind int

const (
	DIT DivisionKind = iota // Division-in-Time
	DIF                     // Division-in-Frequency
)

func log2(n uint) (log uint) {
	for ; n > 1; log, n = log+1, n>>1 {
	}
	return
}
