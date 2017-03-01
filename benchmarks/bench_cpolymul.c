#include "mpfrcx.h"
#include "flint/profiler.h"
#include "arb_poly.h"
#include "acb_poly.h"
#include "mpc.h"

int main()
{
    mpcx_t A, B, C;

    long size, bits;
    long i;

    mpc_t t;

    acb_poly_t X, Y, Z;
    acb_t u;
    acb_init(u);

    for (size = 1000000; size <= 1000000; size *= 10)
    {
        for (bits = 100; bits <= 1000; bits *= 10)
        {
            if (_acb_vec_estimate_allocated_bytes(size, bits) > 3e9)
                continue;

            printf("len=%ld  bits=%ld\n", size, bits);

            mpcx_init(A, size, bits);
            mpcx_init(B, size, bits);
            mpcx_init(C, size, bits);
            mpc_init2(t, bits);

            for (i = 0; i < size; i++)
            {
                mpfr_set_ui(mpc_realref(t), 1, MPFR_RNDD);
                mpfr_div_ui(mpc_realref(t), mpc_realref(t), i + 1, MPFR_RNDD);
                mpfr_set_ui(mpc_imagref(t), 1, MPFR_RNDD);
                mpfr_div_ui(mpc_imagref(t), mpc_imagref(t), i + 2, MPFR_RNDD);
                mpcx_set_coeff(A, i, t);

                mpfr_set_ui(mpc_realref(t), 1, MPFR_RNDD);
                mpfr_div_ui(mpc_realref(t), mpc_realref(t), i + 3, MPFR_RNDD);
                mpfr_set_ui(mpc_imagref(t), 1, MPFR_RNDD);
                mpfr_div_ui(mpc_imagref(t), mpc_imagref(t), i + 4, MPFR_RNDD);
                mpcx_set_coeff(B, i, t);
            }

            TIMEIT_START
            mpcx_mul(C, A, B);
            TIMEIT_STOP

            mpc_clear(t);

            mpcx_clear(A);
            mpcx_clear(B);
            mpcx_clear(C);

            acb_poly_init(X);
            acb_poly_init(Y);
            acb_poly_init(Z);

            acb_poly_fit_length(X, size);
            acb_poly_fit_length(Y, size);

            for (i = 0; i < size; i++)
            {
                arb_set_ui(acb_realref(X->coeffs + i), 1);
                arb_div_ui(acb_realref(X->coeffs + i), acb_realref(X->coeffs + i), i + 1, bits);
                arb_set_ui(acb_imagref(X->coeffs + i), 1);
                arb_div_ui(acb_imagref(X->coeffs + i), acb_imagref(X->coeffs + i), i + 2, bits);

                arb_set_ui(acb_realref(Y->coeffs + i), 1);
                arb_div_ui(acb_realref(Y->coeffs + i), acb_realref(Y->coeffs + i), i + 3, bits);
                arb_set_ui(acb_imagref(Y->coeffs + i), 1);
                arb_div_ui(acb_imagref(Y->coeffs + i), acb_imagref(Y->coeffs + i), i + 4, bits);
            }

            _acb_poly_set_length(X, size);
            _acb_poly_set_length(Y, size);

            TIMEIT_START
            acb_poly_mul(Z, X, Y, bits);
            TIMEIT_STOP
            printf("\n");

            acb_poly_clear(X);
            acb_poly_clear(Y);
            acb_poly_clear(Z);

        }
    }
}

