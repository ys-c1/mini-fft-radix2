#include "fft.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define TEST_PI 3.141592653589793238462643383279502884

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(condition, ...)                    \
    do {                                         \
        ++g_checks;                              \
        if (!(condition)) {                      \
            ++g_failures;                        \
            printf("FAIL check %d: ", g_checks); \
            printf(__VA_ARGS__);                 \
            printf("\n");                        \
        }                                        \
    } while (0)

static uint64_t g_rng_state = UINT64_C(0x243f6a8885a308d3);

static double next_unit(void)
{
    g_rng_state = g_rng_state * UINT64_C(6364136223846793005)
                  + UINT64_C(1442695040888963407);
    return 2.0 * (double)((g_rng_state >> 11) & UINT64_C(0x1fffffffffffff))
               / 9007199254740992.0
           - 1.0;
}

static void test_fft_matches_naive(void)
{
    const int n = 8;
    double src_re[8];
    double src_im[8];
    double fft_re[8];
    double fft_im[8];
    double dft_re[8];
    double dft_im[8];
    int i;

    for (i = 0; i < n; ++i) {
        src_re[i] = next_unit();
        src_im[i] = next_unit();
        fft_re[i] = src_re[i];
        fft_im[i] = src_im[i];
    }

    CHECK(fft_forward(fft_re, fft_im, n) == 0,
          "n=8 fft_forward returned non-zero");
    CHECK(dft_naive(src_re, src_im, dft_re, dft_im, n) == 0,
          "n=8 dft_naive returned non-zero");

    for (i = 0; i < n; ++i) {
        CHECK(fabs(fft_re[i] - dft_re[i]) < 1e-9,
              "n=8 real mismatch at bin %d: fft=%.17g dft=%.17g",
              i, fft_re[i], dft_re[i]);
        CHECK(fabs(fft_im[i] - dft_im[i]) < 1e-9,
              "n=8 imag mismatch at bin %d: fft=%.17g dft=%.17g",
              i, fft_im[i], dft_im[i]);
    }
}

static void test_single_tone(void)
{
    const int n = 64;
    const int expected_bin = 5;
    double re[64];
    double im[64];
    double mag[64];
    int peak = 0;
    int k;

    for (k = 0; k < n; ++k) {
        const double phase =
            2.0 * TEST_PI * (double)expected_bin * (double)k / (double)n;

        re[k] = cos(phase);
        im[k] = sin(phase);
    }

    CHECK(fft_forward(re, im, n) == 0,
          "n=64 single-tone fft_forward returned non-zero");

    for (k = 0; k < n; ++k) {
        mag[k] = sqrt(re[k] * re[k] + im[k] * im[k]);
        if (mag[k] > mag[peak]) {
            peak = k;
        }
    }

    CHECK(peak == expected_bin,
          "single-tone peak at bin %d, expected bin %d",
          peak, expected_bin);
    CHECK(fabs(mag[expected_bin] - (double)n) < 1e-9,
          "single-tone peak magnitude %.17g, expected %d",
          mag[expected_bin], n);

    for (k = 0; k < n; ++k) {
        if (k != expected_bin) {
            CHECK(mag[k] < 1e-6,
                  "single-tone leakage at bin %d: magnitude %.17g",
                  k, mag[k]);
        }
    }
}

static void test_round_trip(void)
{
    const int n = 64;
    double original_re[64];
    double original_im[64];
    double re[64];
    double im[64];
    double max_re_error = 0.0;
    double max_im_error = 0.0;
    int i;

    for (i = 0; i < n; ++i) {
        original_re[i] = next_unit();
        original_im[i] = next_unit();
        re[i] = original_re[i];
        im[i] = original_im[i];
    }

    CHECK(fft_forward(re, im, n) == 0,
          "n=64 round-trip fft_forward returned non-zero");
    CHECK(fft_inverse(re, im, n) == 0,
          "n=64 round-trip fft_inverse returned non-zero");

    for (i = 0; i < n; ++i) {
        const double re_error = fabs(re[i] - original_re[i]);
        const double im_error = fabs(im[i] - original_im[i]);

        if (re_error > max_re_error) {
            max_re_error = re_error;
        }
        if (im_error > max_im_error) {
            max_im_error = im_error;
        }
    }

    CHECK(max_re_error < 1e-12,
          "round-trip max real error %.17g exceeds 1e-12", max_re_error);
    CHECK(max_im_error < 1e-12,
          "round-trip max imag error %.17g exceeds 1e-12", max_im_error);
}

static void test_invalid_lengths(void)
{
    double re[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    double im[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    CHECK(fft_forward(re, im, 3) == -1,
          "fft_forward accepted n=3");
    CHECK(fft_forward(re, im, 6) == -1,
          "fft_forward accepted n=6");
    CHECK(fft_forward(re, im, 0) == -1,
          "fft_forward accepted n=0");
    CHECK(fft_forward(re, im, -8) == -1,
          "fft_forward accepted n=-8");

    CHECK(fft_inverse(re, im, 3) == -1,
          "fft_inverse accepted n=3");
    CHECK(fft_inverse(re, im, 6) == -1,
          "fft_inverse accepted n=6");
    CHECK(fft_inverse(re, im, 0) == -1,
          "fft_inverse accepted n=0");
    CHECK(fft_inverse(re, im, -8) == -1,
          "fft_inverse accepted n=-8");
}

static void test_dc_signal(void)
{
    const int n = 8;
    double re[8];
    double im[8];
    int k;

    for (k = 0; k < n; ++k) {
        re[k] = 1.0;
        im[k] = 0.0;
    }

    CHECK(fft_forward(re, im, n) == 0,
          "n=8 DC fft_forward returned non-zero");
    CHECK(fabs(re[0] - (double)n) < 1e-9,
          "DC bin 0 real part %.17g, expected %d", re[0], n);
    CHECK(fabs(im[0]) < 1e-9,
          "DC bin 0 imag part %.17g, expected 0", im[0]);

    for (k = 1; k < n; ++k) {
        CHECK(fabs(re[k]) < 1e-9,
              "DC bin %d real part %.17g, expected 0", k, re[k]);
    }
    for (k = 1; k < n; ++k) {
        CHECK(fabs(im[k]) < 1e-9,
              "DC bin %d imag part %.17g, expected 0", k, im[k]);
    }
}

int main(void)
{
    printf("mini-fft-radix2 test suite\n");

    test_fft_matches_naive();
    test_single_tone();
    test_round_trip();
    test_invalid_lengths();
    test_dc_signal();

    printf("RESULT: checks=%d failures=%d\n", g_checks, g_failures);

    return g_failures == 0 ? 0 : 1;
}
