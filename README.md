# mini-fft-radix2

A small iterative radix-2 FFT written in C11. It has no third-party
dependencies and uses only the C standard library plus `libm`.

## Files

- `src/fft.h` - public API.
- `src/fft.c` - FFT and naive DFT implementations.
- `tests/test_fft.c` - self-checking test program.
- `Makefile` - build and test entry points.

## API

```c
int fft_forward(double *re, double *im, int n);
int fft_inverse(double *re, double *im, int n);
int dft_naive(const double *re, const double *im,
              double *ore, double *oim, int n);
```

`fft_forward` and `fft_inverse` transform the complex data in place. `n` must
be a positive power of two. `fft_inverse` normalizes its result by `1/n`.
`dft_naive` is an O(n^2) reference DFT and writes to separate output arrays.

All three functions return `0` on success. They return `-1` for invalid
arguments, and the FFT functions also return `-1` if temporary table
allocation fails.

## Build and test

```sh
make test
```

The test runner prints a final line of the form:

```text
RESULT: checks=N failures=M
```

It exits with status 0 when every check passes and status 1 otherwise. Use
`make clean` to remove the build directory.

## Implementation notes

The transform performs an in-place bit-reversal permutation followed by
iterative Cooley-Tukey butterflies. Twiddle factors are precomputed in a
table of `n/2` complex values. The forward transform uses the
`exp(-2*pi*i*k/n)` convention, the inverse uses the conjugate convention,
and the inverse result is scaled by `1/n`.
