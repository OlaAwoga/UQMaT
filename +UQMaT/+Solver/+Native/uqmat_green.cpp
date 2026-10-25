/*
 * uqmat_green.cpp — MEX Green's function solver for UQMaT
 *
 * Computes G(ω) = (ω+iη - H)^{-1} for a batch of ω values using
 * direct matrix inversion via LAPACK zgesv (LU factorisation).
 *
 * Advantages over MATLAB loop:
 *   - Single MEX call overhead instead of Nω calls
 *   - LAPACK zgesv is cache-optimised LU
 *   - Can be extended with OpenMP over ω points
 *
 * Usage from MATLAB (via NativeGreen.m):
 *   [D, L] = uqmat_green_mex(H, omega, eta)
 *   [D, L] = uqmat_green_mex(H, omega, eta, sites)
 *
 * Inputs:
 *   H     : N×N real or complex Hamiltonian (double)
 *   omega : Nw×1 energy grid
 *   eta   : scalar broadening
 *   sites : (optional) 1×Ns site indices for LDOS [1-based]
 *
 * Outputs:
 *   D    : Nw×1 DOS = -1/π Im Tr G(ω)
 *   L    : Nw×Ns LDOS = -1/π Im G(i,i,ω) for sites i
 *
 * Citation:
 *   Anderson E. et al., LAPACK Users' Guide, SIAM (1999).
 *   Economou E.N., Green's Functions in Quantum Physics, Springer (2006).
 */

#include "mex.h"
#include "lapack.h"
#include <cmath>
#include <vector>
#include <complex>

typedef std::complex<double> cdbl;

/* Solve (zI - H) X = I  for complex z, return diagonal of X */
static void green_diag(const std::vector<cdbl>& H, mwSignedIndex N,
                       cdbl z, const std::vector<mwSignedIndex>& sites,
                       double& dos_out, std::vector<double>& ldos_out)
{
    mwSignedIndex N2 = N * N;
    /* Build A = zI - H as interleaved double array for LAPACK */
    std::vector<double> A(2 * N2);
    for (mwSignedIndex j = 0; j < N; ++j)
        for (mwSignedIndex i = 0; i < N; ++i) {
            cdbl val = (i == j ? z : cdbl(0,0)) - H[j*N+i];
            A[2*(j*N+i)]   = val.real();
            A[2*(j*N+i)+1] = val.imag();
        }

    /* RHS = identity */
    std::vector<double> B(2 * N2, 0.0);
    for (mwSignedIndex i = 0; i < N; ++i) {
        B[2*(i*N+i)]   = 1.0;
        B[2*(i*N+i)+1] = 0.0;
    }

    std::vector<mwSignedIndex> ipiv(N);
    mwSignedIndex info = 0;

    zgesv(&N, &N, A.data(), &N, ipiv.data(), B.data(), &N, &info);
    if (info != 0) { dos_out = 0; return; }

    /* DOS = -1/π Im Tr G */
    double tr_im = 0.0;
    for (mwSignedIndex i = 0; i < N; ++i)
        tr_im += B[2*(i*N+i)+1];
    dos_out = -tr_im / (M_PI * N);

    /* LDOS at requested sites */
    for (mwSize si = 0; si < sites.size(); ++si) {
        mwSignedIndex idx = sites[si] - 1;  /* 0-based */
        if (idx >= 0 && idx < N)
            ldos_out[si] = -B[2*(idx*N+idx)+1] / M_PI;
    }
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (nrhs < 3)
        mexErrMsgIdAndTxt("UQMaT:green:args",
            "Usage: [D,L] = uqmat_green_mex(H, omega, eta [, sites])");

    /* Parse inputs */
    mwSignedIndex N = (mwSignedIndex)mxGetM(prhs[0]);
    double* Hr = mxGetPr(prhs[0]);
    double* Hi = mxIsComplex(prhs[0]) ? mxGetPi(prhs[0]) : nullptr;

    double* omega_ptr = mxGetPr(prhs[1]);
    mwSize  Nw        = mxGetNumberOfElements(prhs[1]);
    double  eta       = mxGetScalar(prhs[2]);

    std::vector<mwSignedIndex> sites;
    if (nrhs >= 4) {
        double* sp = mxGetPr(prhs[3]);
        mwSize  Ns = mxGetNumberOfElements(prhs[3]);
        sites.resize(Ns);
        for (mwSize k = 0; k < Ns; ++k)
            sites[k] = (mwSignedIndex)sp[k];
    }

    /* Build complex H */
    std::vector<cdbl> H(N*N, cdbl(0,0));
    for (mwSignedIndex k = 0; k < N*N; ++k)
        H[k] = cdbl(Hr[k], Hi ? Hi[k] : 0.0);

    /* Outputs */
    plhs[0] = mxCreateDoubleMatrix(Nw, 1, mxREAL);
    double* D = mxGetPr(plhs[0]);

    bool do_ldos = (nlhs > 1 && !sites.empty());
    mxArray* mL = nullptr;
    double*  Lp = nullptr;
    if (do_ldos) {
        mL = mxCreateDoubleMatrix(Nw, (mwSize)sites.size(), mxREAL);
        Lp = mxGetPr(mL);
    }

    /* Loop over ω */
    for (mwSize iw = 0; iw < Nw; ++iw) {
        cdbl z(omega_ptr[iw], eta);
        std::vector<double> ldos_row(sites.size(), 0.0);
        green_diag(H, N, z, sites, D[iw], ldos_row);
        if (do_ldos)
            for (mwSize si = 0; si < sites.size(); ++si)
                Lp[iw + si*Nw] = ldos_row[si];
    }

    if (do_ldos) plhs[1] = mL;
}
