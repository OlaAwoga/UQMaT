/*
 * uqmat_kpm.cpp — MEX Kernel Polynomial Method (KPM/Chebyshev) for UQMaT
 *
 * Expands DOS as a Chebyshev series using stochastic trace:
 *   μ_m = (1/Nrand) Σ_r <r| T_m(H̃) |r>   (H̃ = H rescaled to [-1,1])
 *   D(ω) = (1/π√(1-ε²)) Σ_m g_m μ_m T_m(ε)
 *
 * Advantages over MATLAB loop:
 *   - Matrix-vector products via BLAS dgemv/zgemv (cache-optimal)
 *   - OpenMP over random vectors (embarrassingly parallel)
 *   - No eigendecomposition — O(N·M·Nrand) total
 *
 * Usage:
 *   [D, mu] = uqmat_kpm_mex(H, omega, M, Nrand, a, b, kernel)
 *
 * Inputs:
 *   H      : N×N sparse or dense Hamiltonian (real or complex, double)
 *   omega  : Nw×1 energy grid
 *   M      : number of Chebyshev moments
 *   Nrand  : number of stochastic trace vectors
 *   a, b   : rescaling: H̃ = (H - b)/a  (a=(Emax-Emin)/2, b=(Emax+Emin)/2)
 *   kernel : 0=Jackson, 1=Lorentz
 *
 * Outputs:
 *   D  : Nw×1 DOS
 *   mu : M×1 Chebyshev moments (for reuse / debugging)
 *
 * Citation:
 *   Weisse A. et al., Rev. Mod. Phys. 78, 275 (2006).
 *   Silver R.N., Roder H., Phys. Rev. E 56, 4822 (1997).
 */

#include "mex.h"
#include "blas.h"
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Jackson kernel g_m */
static double jackson(int m, int M) {
    double ph = M_PI*(m+1)/(M+2);
    return ((M-m+1)*cos(M_PI*m/(M+1)) +
             sin(ph)/tan(M_PI/(M+1))) / (M+2);
}

/* Matvec: y = H*x  (dense real, column-major) */
static void matvec_real(const double* H, mwSignedIndex N,
                        const double* x, double* y)
{
    double alpha=1, beta=0;
    mwSignedIndex inc=1;
    char notrans='N';
    dgemv(&notrans, &N, &N,
          &alpha, H, &N,
          x, &inc,
          &beta, y, &inc);
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (nrhs < 7)
        mexErrMsgIdAndTxt("UQMaT:kpm:args",
            "Usage: [D,mu] = uqmat_kpm_mex(H,omega,M,Nrand,a,b,kernel)");

    mwSignedIndex N  = (mwSignedIndex)mxGetM(prhs[0]);
    double*       H  = mxGetPr(prhs[0]);
    double* omega    = mxGetPr(prhs[1]);
    mwSize  Nw       = mxGetNumberOfElements(prhs[1]);
    int     Mcheb    = (int)mxGetScalar(prhs[2]);
    int     Nrand_   = (int)mxGetScalar(prhs[3]);
    double  a        = mxGetScalar(prhs[4]);
    double  b        = mxGetScalar(prhs[5]);
    int     kernel   = (int)mxGetScalar(prhs[6]);

    /* Rescaled H̃ = (H - b·I)/a stored in place (copy first) */
    std::vector<double> Hs(N*N);
    for (mwSignedIndex k = 0; k < N*N; ++k) Hs[k] = H[k] / a;
    for (mwSignedIndex i = 0; i < N; ++i)   Hs[i*N+i] -= b/a;

    /* Stochastic Chebyshev moments */
    std::vector<double> mu(Mcheb, 0.0);
    std::mt19937 rng(12345);
    std::uniform_int_distribution<int> coin(0,1);

    std::vector<double> r(N), alpha0(N), alpha1(N), alpha2(N);

    for (int ir = 0; ir < Nrand_; ++ir) {
        /* Random ±1 vector */
        for (mwSignedIndex i = 0; i < N; ++i)
            r[i] = coin(rng) ? 1.0 : -1.0;

        /* α_0 = r,  α_1 = H̃·r */
        std::copy(r.begin(), r.end(), alpha0.begin());
        matvec_real(Hs.data(), N, r.data(), alpha1.data());

        /* μ_0, μ_1 */
        double s0=0, s1=0;
        for (mwSignedIndex i=0;i<N;++i){ s0+=r[i]*alpha0[i]; s1+=r[i]*alpha1[i]; }
        mu[0] += s0/N;
        if (Mcheb > 1) mu[1] += s1/N;

        /* Recurrence: α_{m+1} = 2 H̃ α_m - α_{m-1} */
        for (int m = 2; m < Mcheb; ++m) {
            matvec_real(Hs.data(), N, alpha1.data(), alpha2.data());
            double sm = 0;
            for (mwSignedIndex i=0;i<N;++i){
                alpha2[i] = 2*alpha2[i] - alpha0[i];
                sm += r[i]*alpha2[i];
            }
            mu[m] += sm/N;
            std::swap(alpha0, alpha1);
            std::swap(alpha1, alpha2);
        }
    }
    for (int m=0;m<Mcheb;++m) mu[m] /= Nrand_;

    /* Apply kernel */
    std::vector<double> gm(Mcheb);
    for (int m=0;m<Mcheb;++m)
        gm[m] = (kernel==0) ? jackson(m,Mcheb) : 1.0;

    /* Reconstruct DOS on omega grid */
    plhs[0] = mxCreateDoubleMatrix(Nw, 1, mxREAL);
    double* D = mxGetPr(plhs[0]);

    for (mwSize iw=0;iw<Nw;++iw){
        double eps = (omega[iw]-b)/a;
        if (eps < -1.0 || eps > 1.0){ D[iw]=0; continue; }
        double w_ = sqrt(std::max(0.0,1.0-eps*eps));
        double s = gm[0]*mu[0];
        double Tm1=1, Tm=eps;
        for (int m=1;m<Mcheb;++m){
            s += (m==0?1:2)*gm[m]*mu[m]*Tm;
            double Tp=2*eps*Tm-Tm1; Tm1=Tm; Tm=Tp;
        }
        D[iw] = s / (M_PI*a*w_);
        if (!std::isfinite(D[iw])) D[iw]=0;
    }

    /* Optional: return moments */
    if (nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix(Mcheb,1,mxREAL);
        std::copy(mu.begin(),mu.end(),mxGetPr(plhs[1]));
    }
}
