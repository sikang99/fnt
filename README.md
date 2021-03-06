# Fast Numeric Transforms for Go
Fast Numeric Transforms for Go. Provides 1-dimensional discrete Hartley and Fourier transforms. Both transforms are implemented as iterative radix-2 Cooley-Tukey algorithms. Work is based loosely on implementations given by Jörg Arndt in [Matters Computational](http://www.jjj.de/fxt/#fxtbook).

## Example

```Go
package main

import (
	"fmt"
	"github.com/bemasher/fnt"
)

func main() {
	fft := fnt.NewFFT(8)

	samples := []complex128{1, 1, 1, 1, 0, 0, 0, 0}
	fmt.Printf("%+0.3f\n", samples)

	// Transform to Fourier domain.
	fft.Execute(samples, fnt.DIT, fnt.Forward, false)
	fmt.Printf("%+0.3f\n", samples)

	// Transform back to time-domain and normalize.
	fft.Execute(samples, fnt.DIT, fnt.Backward, true)
	fmt.Printf("%+0.3f\n", samples)
}
```

Output:
```
[(+1.000+0.000i) (+1.000+0.000i) (+1.000+0.000i) (+1.000+0.000i) (+0.000+0.000i) (+0.000+0.000i) (+0.000+0.000i) (+0.000+0.000i)]
[(+4.000+0.000i) (+1.000-2.414i) (+0.000+0.000i) (+1.000-0.414i) (+0.000+0.000i) (+1.000+0.414i) (+0.000+0.000i) (+1.000+2.414i)]
[(+1.000+0.000i) (+1.000+0.000i) (+1.000-0.000i) (+1.000-0.000i) (+0.000+0.000i) (+0.000-0.000i) (+0.000+0.000i) (+0.000+0.000i)]
```

## Documentation

[![GoDoc](https://godoc.org/github.com/bemasher/fnt?status.png)](https://godoc.org/github.com/bemasher/fnt)

## To Do
 * Explore more localized algorithm for both FFT and FHT.
 * Explore parallelization for all transforms.
 * Implement convolution methods for both FFT and FHT.
 * Implement Real-to-Complex and Complex-to-Real FFT's.
 * Implement a Number Theoretic Transform for exact convolution.
 * Implement tests using relationships between FFT and FHT.
