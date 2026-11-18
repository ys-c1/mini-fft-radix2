#ifndef MINI_FFT_RADIX2_FFT_H
#define MINI_FFT_RADIX2_FFT_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * In-place radix-2 FFT.
 *
 * re/im hold the complex input and receive the complex output. n must be a
 * power of two and both pointers must be non-NULL.
 *
 * Returns 0 on success, -1 for invalid arguments or allocation failure.
 */
int fft_forward(double *re, double *im, int n);

/*
 * In-place inverse radix-2 FFT. The result is normalized by 1/n.
 *
 * Returns 0 on success, -1 for invalid arguments or allocation failure.
 */
int fft_inverse(double *re, double *im, int n);

/*
 * Naive O(n^2) DFT used as a reference implementation.
 *
 * ore/oim must not alias re/im. Returns 0 on success, -1 for invalid
 * arguments.
 */
int dft_naive(const double *re, const double *im,
              double *ore, double *oim, int n);

#ifdef __cplusplus
}
#endif

#endif /* MINI_FFT_RADIX2_FFT_H */
