#include "fft.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define FFT_PI 3.141592653589793238462643383279502884

static int is_power_of_two(int n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

static unsigned int reverse_bits(unsigned int value, int bits)
{
    unsigned int result = 0u;
    int i;

    for (i = 0; i < bits; ++i) {
        result = (result << 1) | (value & 1u);
        value >>= 1;
    }

    return result;
}

static void bit_reverse_permute(double *re, double *im, int n)
{
    int bits = 0;
    int i;

    for (i = n; i > 1; i >>= 1) {
        ++bits;
    }

    for (i = 0; i < n; ++i) {
        unsigned int j = reverse_bits((unsigned int)i, bits);

        if ((int)j > i) {
            double tr = re[i];
            double ti = im[i];

            re[i] = re[(int)j];
            im[i] = im[(int)j];
            re[(int)j] = tr;
            im[(int)j] = ti;
        }
    }
}

static int build_twiddle_table(int n, double **tw_re, double **tw_im)
{
    const int half = n / 2;
    double *table_re;
    double *table_im;
    int k;

    *tw_re = NULL;
    *tw_im = NULL;

    if (half == 0) {
        return 0;
    }

    table_re = (double *)malloc((size_t)half * sizeof(*table_re));
    table_im = (double *)malloc((size_t)half * sizeof(*table_im));
    if (table_re == NULL || table_im == NULL) {
        free(table_re);
        free(table_im);
        return -1;
    }

    for (k = 0; k < half; ++k) {
        const double angle = 2.0 * FFT_PI * (double)k / (double)n;

        table_re[k] = cos(angle);
        table_im[k] = sin(angle);
    }

    *tw_re = table_re;
    *tw_im = table_im;
    return 0;
}

static void fft_transform(double *re, double *im, int n,
                          const double *tw_re, const double *tw_im,
                          int inverse)
{
    int len;

    bit_reverse_permute(re, im, n);

    for (len = 2; len <= n; len <<= 1) {
        const int half = len >> 1;
        const int stride = n / len;
        int start;

        for (start = 0; start < n; start += len) {
            int j;

            for (j = 0; j < half; ++j) {
                const int t = j * stride;
                const int a = start + j;
                const int b = a + half;
                const double wr = tw_re[t];
                const double wi = inverse ? tw_im[t] : tw_im[t];
                const double ur = re[a];
                const double ui = im[a];
                const double vr = re[b] * wr - im[b] * wi;
                const double vi = re[b] * wi + im[b] * wr;

                re[a] = ur + vr;
                im[a] = ui + vi;
                re[b] = ur - vr;
                im[b] = ui - vi;
            }
        }
    }
}

int fft_forward(double *re, double *im, int n)
{
    double *tw_re;
    double *tw_im;

    if (re == NULL || im == NULL || !is_power_of_two(n)) {
        return -1;
    }

    if (n == 1) {
        return 0;
    }

    if (build_twiddle_table(n, &tw_re, &tw_im) != 0) {
        return -1;
    }

    fft_transform(re, im, n, tw_re, tw_im, 0);

    free(tw_re);
    free(tw_im);
    return 0;
}

int fft_inverse(double *re, double *im, int n)
{
    double *tw_re;
    double *tw_im;

    if (re == NULL || im == NULL || !is_power_of_two(n)) {
        return -1;
    }

    if (n == 1) {
        return 0;
    }

    if (build_twiddle_table(n, &tw_re, &tw_im) != 0) {
        return -1;
    }

    fft_transform(re, im, n, tw_re, tw_im, 1);


    free(tw_re);
    free(tw_im);

    return 0;
}

int dft_naive(const double *re, const double *im,
              double *ore, double *oim, int n)
{
    int k;

    if (re == NULL || im == NULL || ore == NULL || oim == NULL || n <= 0) {
        return -1;
    }

    for (k = 0; k < n; ++k) {
        double sum_re = 0.0;
        double sum_im = 0.0;
        int t;

        for (t = 0; t < n; ++t) {
            const double angle =
                -2.0 * FFT_PI * (double)k * (double)t / (double)n;
            const double c = cos(angle);
            const double s = sin(angle);

            sum_re += re[t] * c - im[t] * s;
            sum_im += re[t] * s + im[t] * c;
        }

        ore[k] = sum_re;
        oim[k] = sum_im;
    }

    return 0;
}
