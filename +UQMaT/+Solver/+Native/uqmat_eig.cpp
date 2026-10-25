/*
 * uqmat_eig.cpp — MEX eigensolver for UQMaT
 *
 * Wraps LAPACK dsyev/zheev (real/complex Hermitian) via MATLAB's built-in
 * LAPACK interface. Faster than MATLAB eig() for large systems because:
 *   - Zero MATLAB overhead in the inner loop
 *   - Can use OpenMP for parallel BLAS (if linked)
 *   - Direct LAPACK call, no intermediate copies
 *
 * Compilation (from MATLAB):
 *   UQMaT.build_mex()     % automatic
 *   % or manually:
 *   mex -O COPTIMFLAGS='-O3 -march=native' uqmat_eig.cpp -llapack -lblas
 *
 * Usage from MATLAB (called automatically via NativeEig.m):
 *   [E, V] = uqmat_eig_mex(H)         % real symmetric H
 *   [E, V] = uqmat_eig_mex(H, 'vec')  % real, with eigenvectors
 *   [E, V] = uqmat_eig_mex(H, 'novec')% real, eigenvalues only (faster)
 *
 * Inputs:
 *   H    : N×N real symmetric or complex Hermitian matrix
 *   mode : 'vec' (default) | 'novec'
 *
 * Outputs:
 *   E    : N×1 sorted eigenvalues (ascending)
 *   V    : N×N eigenvectors (columns), or [] if 'novec'
 *
 * Citation:
 *   Anderson E. et al., LAPACK Users' Guide, SIAM (1999).
 *   https://www.netlib.org/lapack/
 */

#include "mex.h"
#include "lapack.h"
#include <cstring>
#include <algorithm>
#include <vector>

/* ── Real symmetric: dsyev ─────────────────────────────────────────────────── */
static void eig_real(double* H, mwSignedIndex N, bool compute_vec,
                     double* E, double* V)
{
    char jobz = compute_vec ? 'V' : 'N';
    char uplo = 'U';
    mwSignedIndex lda = N;
    mwSignedIndex lwork = -1, info = 0;
    double wkopt = 0;

    /* Copy H into V (dsyev overwrites) */
    if (compute_vec)
        std::copy(H, H + N*N, V);
    else
        std::copy(H, H + N*N, E);  /* temp use of E buffer */

    double* A = compute_vec ? V : E;

    /* Query workspace */
    dsyev(&jobz, &uplo, &N, A, &lda, E, &wkopt, &lwork, &info);
    lwork = (mwSignedIndex)wkopt;
    std::vector<double> work(lwork);

    /* Compute */
    dsyev(&jobz, &uplo, &N, A, &lda, E, work.data(), &lwork, &info);

    if (info != 0)
        mexErrMsgIdAndTxt("UQMaT:eig:lapack",
            "dsyev failed with info=%d", (int)info);
}

/* ── Complex Hermitian: zheev ──────────────────────────────────────────────── */
static void eig_complex(double* Hr, double* Hi, mwSignedIndex N, bool compute_vec,
                        double* E, double* Vr, double* Vi)
{
    char jobz = compute_vec ? 'V' : 'N';
    char uplo = 'U';
    mwSignedIndex lda = N;
    mwSignedIndex lwork = -1, info = 0;
    double rwork_query[1];
    mwSignedIndex lrwork = 3*N - 2;
    std::vector<double> rwork(lrwork);

    /* Interleave real/imag into complex array */
    mwSize n2 = (mwSize)(N * N);
    std::vector<double> A(2*n2);
    for (mwSize k = 0; k < n2; ++k) {
        A[2*k]   = Hr[k];
        A[2*k+1] = Hi[k];
    }

    /* Query workspace */
    double wkopt[2] = {0,0};
    zheev(&jobz, &uplo, &N, A.data(), &lda, E,
          wkopt, &lwork, rwork.data(), &info);
    lwork = (mwSignedIndex)wkopt[0];
    std::vector<double> work(2*lwork);

    /* Compute */
    zheev(&jobz, &uplo, &N, A.data(), &lda, E,
          work.data(), &lwork, rwork.data(), &info);

    if (info != 0)
        mexErrMsgIdAndTxt("UQMaT:eig:lapack",
            "zheev failed with info=%d", (int)info);

    /* Deinterleave back into real/imag */
    if (compute_vec && Vr && Vi) {
        for (mwSize k = 0; k < n2; ++k) {
            Vr[k] = A[2*k];
            Vi[k] = A[2*k+1];
        }
    }
}

/* ── MEX gateway ───────────────────────────────────────────────────────────── */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (nrhs < 1)
        mexErrMsgIdAndTxt("UQMaT:eig:args","Usage: [E,V] = uqmat_eig_mex(H)");

    const mxArray* mH = prhs[0];
    if (!mxIsDouble(mH) || mxGetNumberOfDimensions(mH) != 2)
        mexErrMsgIdAndTxt("UQMaT:eig:type","H must be a double matrix.");

    mwSignedIndex N = (mwSignedIndex)mxGetM(mH);
    if (N != (mwSignedIndex)mxGetN(mH))
        mexErrMsgIdAndTxt("UQMaT:eig:square","H must be square.");

    bool compute_vec = (nlhs > 1);
    bool is_complex  = mxIsComplex(mH);

    /* Output: E (N×1) */
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double* E = mxGetPr(plhs[0]);

    if (!is_complex) {
        double* Hr = mxGetPr(mH);
        if (compute_vec) {
            plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
            eig_real(Hr, N, true, E, mxGetPr(plhs[1]));
        } else {
            eig_real(Hr, N, false, E, nullptr);
        }
    } else {
        double* Hr = mxGetPr(mH);
        double* Hi = mxGetPi(mH);
        if (compute_vec) {
            plhs[1] = mxCreateDoubleMatrix(N, N, mxCOMPLEX);
            eig_complex(Hr, Hi, N, true, E,
                        mxGetPr(plhs[1]), mxGetPi(plhs[1]));
        } else {
            eig_complex(Hr, Hi, N, false, E, nullptr, nullptr);
        }
    }
}
