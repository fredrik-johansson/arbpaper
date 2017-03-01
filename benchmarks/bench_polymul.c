#include "mpfrcx.h"
#include "flint/profiler.h"
#include "arb_poly.h"

int main()
{
    mpfrx_t A, B, C;

    long size, bits;
    long i;

    mpfr_t t;

    arb_poly_t X, Y, Z;
    arb_t u;

    arb_init(u);

    for (size = 1000000; size <= 1000000; size *= 10)
    {
        for (bits = 100; bits <= 100000; bits *= 10)
        {
            if (_arb_vec_estimate_allocated_bytes(size, bits) > 3e9)
                continue;

            printf("len=%ld  bits=%ld    ", size, bits);
            printf("\n");

            mpfrx_init(A, size, bits);
            mpfrx_init(B, size, bits);
            mpfrx_init(C, size, bits);

            mpfr_init2(t, bits);

            for (i = 0; i < size; i++)
            {
                mpfr_set_ui(t, 1, MPFR_RNDD);
                mpfr_div_ui(t, t, i + 1, MPFR_RNDD);
                mpfrx_set_coeff(A, i, t);
                mpfr_set_ui(t, 1, MPFR_RNDD);
                mpfr_div_ui(t, t, i + 2, MPFR_RNDD);
                mpfrx_set_coeff(B, i, t);
            }

            TIMEIT_START
            mpfrx_mul(C, A, B);
            TIMEIT_STOP

            mpfr_clear(t);
            mpfrx_clear(A);
            mpfrx_clear(B);
            mpfrx_clear(C);

            arb_poly_init(X);
            arb_poly_init(Y);
            arb_poly_init(Z);

            arb_poly_fit_length(X, size);
            arb_poly_fit_length(Y, size);

            _arb_poly_set_length(X, size);
            _arb_poly_set_length(Y, size);

            for (i = 0; i < size; i++)
            {
                arb_set_ui(X->coeffs + i, 1);
                arb_div_ui(X->coeffs + i, X->coeffs + i, i + 1, bits);
                arb_set_ui(Y->coeffs + i, 1);
                arb_div_ui(Y->coeffs + i, Y->coeffs + i, i + 2, bits);
            }

            TIMEIT_START
            arb_poly_mul(Z, X, Y, bits);
            TIMEIT_STOP

            arb_poly_clear(X);
            arb_poly_clear(Y);
            arb_poly_clear(Z);
        }
    }
}

