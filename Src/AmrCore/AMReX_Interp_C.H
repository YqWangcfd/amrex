#ifndef AMREX_INTERP_C_H_
#define AMREX_INTERP_C_H_
#include <AMReX_Config.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_Interp_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_Interp_2D_C.H>
#else
#include <AMReX_Interp_3D_C.H>
#endif
// #include <AMReX_Print.H>

namespace amrex {

//
// Fill fine values with piecewise-constant interpolation of coarse data.
// Operate only on faces that overlap -- i.e., only fill the fine faces that
// make up each coarse face, leave the in-between faces alone.
//
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_x (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int ci = amrex::coarsen(fi, ratio[0]);
    if (ci*ratio[0] == fi) {
#if (AMREX_SPACEDIM >= 2)
        int cj = amrex::coarsen(fj, ratio[1]);
#else
        int cj = 0;
#endif
#if (AMREX_SPACEDIM == 3)
        int ck = amrex::coarsen(fk, ratio[2]);
#else
        int ck = 0;
#endif

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_y (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int cj = amrex::coarsen(fj, ratio[1]);
    if (cj*ratio[1] == fj) {
        int ci = amrex::coarsen(fi, ratio[0]);
#if (AMREX_SPACEDIM == 3)
        int ck = amrex::coarsen(fk, ratio[2]);
#else
        int ck = 0;
#endif

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_z (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int ck = amrex::coarsen(fk, ratio[2]);
    if (ck*ratio[2] == fk) {
        int ci = amrex::coarsen(fi, ratio[0]);
        int cj = amrex::coarsen(fj, ratio[1]);

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

//
// Fill fine values with tangential interpolation of coarse data.
// Operate only on faces that overlap -- i.e., only fill the fine faces that
// make up each coarse face, leave the in-between faces alone.
//
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_cons_linear_face_interp (int i, int j, int k, int n, Array4<T> const& fine,
                              Array4<T const> const& crse, Array4<int const> const& mask,
                              IntVect const& ratio,  Box const& per_grown_domain, int dim) noexcept
{
    int ci = amrex::coarsen(i, ratio[0]);

#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(per_grown_domain);
    int cj = 0;
#else
    int cj = amrex::coarsen(j, ratio[1]);
#endif

#if (AMREX_SPACEDIM == 3)
    int ck = amrex::coarsen(k, ratio[2]);
#else
    int ck = 0;
#endif

    if (dim == 0 && ci*ratio[0] == i) {
        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(i, j, k, n) = crse(ci, cj, ck, n);
#if (AMREX_SPACEDIM >= 2)
            if (cj > per_grown_domain.smallEnd(1) && cj < per_grown_domain.bigEnd(1) && ratio[1] > 1) {
                Real sfy = Real(1.0);
                Real dc = Real(0.5) * (crse(ci,cj+1,ck,n) - crse(ci,cj-1,ck,n));
                Real df = Real(2.0) * (crse(ci,cj+1,ck,n) - crse(ci,cj  ,ck,n));
                Real db = Real(2.0) * (crse(ci,cj  ,ck,n) - crse(ci,cj-1,ck,n));
                Real sy = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfy = amrex::min(sfy, sy / dc);
                }
                Real slope = dc;
                Real yoff = (static_cast<Real>(j - cj*ratio[1]) + Real(0.5)) / Real(ratio[1]) - Real(0.5);
                fine(i,j,k,n) += yoff * slope * sfy;
            } // jc
#if (AMREX_SPACEDIM == 3)
            if (ck > per_grown_domain.smallEnd(2) && ck < per_grown_domain.bigEnd(2) && ratio[2] > 1) {
                Real sfz = Real(1.0);
                Real dc = Real(0.5) * (crse(ci,cj,ck+1,n) - crse(ci,cj,ck-1,n));
                Real df = Real(2.0) * (crse(ci,cj,ck+1,n) - crse(ci,cj,ck  ,n));
                Real db = Real(2.0) * (crse(ci,cj,ck  ,n) - crse(ci,cj,ck-1,n));
                Real sz = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sz = std::copysign(Real(1.),dc)*amrex::min(sz,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfz = amrex::min(sfz, sz / dc);
                }
                Real slope = dc;
                Real zoff = (static_cast<Real>(k - ck*ratio[2]) + Real(0.5)) / Real(ratio[2]) - Real(0.5);
                fine(i,j,k,n) += zoff * slope * sfz;
            } // ck
#endif
#endif
        } // mask
    } // dim

#if (AMREX_SPACEDIM >= 2)
    if (dim == 1 && cj*ratio[1] == j) {
        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(i, j, k, n) = crse(ci, cj, ck, n);
            if (ci > per_grown_domain.smallEnd(0) && ci < per_grown_domain.bigEnd(0) && ratio[0] > 1) {
                Real sfx = Real(1.0);
                Real dc = Real(0.5) * (crse(ci+1,cj,ck,n) - crse(ci-1,cj,ck,n));
                Real df = Real(2.0) * (crse(ci+1,cj,ck,n) - crse(ci  ,cj,ck,n));
                Real db = Real(2.0) * (crse(ci  ,cj,ck,n) - crse(ci-1,cj,ck,n));
                Real sx = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfx = amrex::min(sfx, sx / dc);
                }
                Real slope = dc;
                Real xoff = (static_cast<Real>(i - ci*ratio[0]) + Real(0.5)) / Real(ratio[0]) - Real(0.5);
                fine(i,j,k,n) += xoff * slope * sfx;
            } // ci
#if (AMREX_SPACEDIM == 3)
            if (ck > per_grown_domain.smallEnd(2) && ck < per_grown_domain.bigEnd(2) && ratio[2] > 1) {
                Real sfz = Real(1.0);
                Real dc = Real(0.5) * (crse(ci,cj,ck+1,n) - crse(ci,cj,ck-1,n));
                Real df = Real(2.0) * (crse(ci,cj,ck+1,n) - crse(ci,cj,ck  ,n));
                Real db = Real(2.0) * (crse(ci,cj,ck  ,n) - crse(ci,cj,ck-1,n));
                Real sz = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sz = std::copysign(Real(1.),dc)*amrex::min(sz,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfz = amrex::min(sfz, sz / dc);
                }
                Real slope = dc;
                Real zoff = (static_cast<Real>(k - ck*ratio[2]) + Real(0.5)) / Real(ratio[2]) - Real(0.5);
                fine(i,j,k,n) += zoff * slope * sfz;
            } // ck
#endif // SPACEDIM >= 3
        } // mask
    } // dim == 1
#endif // SPACEDIM >= 2

#if (AMREX_SPACEDIM == 3)
    if (dim == 2 && ck*ratio[2] == k) {
        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(i, j, k, n) = crse(ci, cj, ck, n);
            if (ci > per_grown_domain.smallEnd(0) && ci < per_grown_domain.bigEnd(0) && ratio[0] > 1) {
                Real sfx = Real(1.0);
                Real dc = Real(0.5) * (crse(ci+1,cj,ck,n) - crse(ci-1,cj,ck,n));
                Real df = Real(2.0) * (crse(ci+1,cj,ck,n) - crse(ci  ,cj,ck,n));
                Real db = Real(2.0) * (crse(ci  ,cj,ck,n) - crse(ci-1,cj,ck,n));
                Real sx = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfx = amrex::min(sfx, sx / dc);
                }
                Real slope = dc;
                Real xoff = (static_cast<Real>(i - ci*ratio[0]) + Real(0.5)) / Real(ratio[0]) - Real(0.5);
                fine(i,j,k,n) += xoff * slope * sfx;
            } // ci
            if (cj > per_grown_domain.smallEnd(1) && cj < per_grown_domain.bigEnd(1) && ratio[1] > 1) {
                Real sfy = Real(1.0);
                Real dc = Real(0.5) * (crse(ci,cj+1,ck,n) - crse(ci,cj-1,ck,n));
                Real df = Real(2.0) * (crse(ci,cj+1,ck,n) - crse(ci,cj  ,ck,n));
                Real db = Real(2.0) * (crse(ci,cj  ,ck,n) - crse(ci,cj-1,ck,n));
                Real sy = (df*db >= Real(0.0)) ?
                    amrex::min(std::abs(df),std::abs(db)) : Real(0.);
                sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));
                if (dc != Real(0.0)) {
                    sfy = amrex::min(sfy, sy / dc);
                }
                Real slope = dc;
                Real yoff = (static_cast<Real>(j - cj*ratio[1]) + Real(0.5)) / Real(ratio[1]) - Real(0.5);
                fine(i,j,k,n) += yoff * slope * sfy;
            } // cj
        } // mask
    } // dim == 2
#endif
}

//
// Do linear in dir, pc transverse to dir, leave alone the fine values
// lining up with coarse edges--assume these have been set to hold the
// values you want to interpolate to the rest.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int ci = amrex::coarsen(i,ratio[0]);

    if (i-ci*ratio[0] != 0)
    {
        Real const w = static_cast<Real>(i-ci*ratio[0]) * (Real(1.)/Real(ratio[0]));
        int i1 = ci*ratio[0];
        int i2 = (ci+1)*ratio[0];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i1,j,k,n) + w * fine(i2,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int cj = amrex::coarsen(j,ratio[1]);

    if (j-cj*ratio[1] != 0)
    {
        Real const w = static_cast<Real>(j-cj*ratio[1]) * (Real(1.)/Real(ratio[1]));
        int j1 = cj*ratio[1];
        int j2 = (cj+1)*ratio[1];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i,j1,k,n) + w * fine(i,j2,k,n);
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_z (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int ck = amrex::coarsen(k,ratio[2]);

    if (k-ck*ratio[2] != 0)
    {
        Real const w = static_cast<Real>(k-ck*ratio[2]) * (Real(1.)/Real(ratio[2]));
        int k1 = ck*ratio[2];
        int k2 = (ck+1)*ratio[2];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i,j,k1,n) + w * fine(i,j,k2,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_x (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984375), Real(-0.123046875),
                                      Real(0.9228515625), Real(0.205078125),
                                      Real(-0.02197265625)};
    int ii = amrex::coarsen(i,2);
    int s = 2*(i-ii*2) - 1;  // if i == ii*2, s = -1; if i == ii*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(ii-2,j,k,n)
        +           c(  -s)*crse(ii-1,j,k,n)
        +           c(   0)*crse(ii  ,j,k,n)
        +           c(   s)*crse(ii+1,j,k,n)
        +           c( 2*s)*crse(ii+2,j,k,n);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_y (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984375), Real(-0.123046875),
                                      Real(0.9228515625), Real(0.205078125),
                                      Real(-0.02197265625)};
    int jj = amrex::coarsen(j,2);
    int s = 2*(j-jj*2) - 1;  // if j == jj*2, s = -1; if j == jj*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(i,jj-2,k,n)
        +           c(  -s)*crse(i,jj-1,k,n)
        +           c(   0)*crse(i,jj  ,k,n)
        +           c(   s)*crse(i,jj+1,k,n)
        +           c( 2*s)*crse(i,jj+2,k,n);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_z (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984375), Real(-0.123046875),
                                      Real(0.9228515625), Real(0.205078125),
                                      Real(-0.02197265625)};
    int kk = amrex::coarsen(k,2);
    int s = 2*(k-kk*2) - 1;  // if k == kk*2, s = -1; if k == kk*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(i,j,kk-2,n)
        +           c(  -s)*crse(i,j,kk-1,n)
        +           c(   0)*crse(i,j,kk  ,n)
        +           c(   s)*crse(i,j,kk+1,n)
        +           c( 2*s)*crse(i,j,kk+2,n);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Lagrange5_interp_x (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    // // 5-point upwind biased
    // constexpr Array1D<Real,-2,2> d = {Real(0.01708984375), Real(-0.107421875),
    //                                   Real(0.3759765625), Real(0.751953125),
    //                                   Real(-0.03759765625)};
    // constexpr Array1D<Real,-2,2> c = {Real(0.01708984375), Real(-0.123046875),
    //                                   Real(0.9228515625), Real(0.205078125),
    //                                   Real(-0.02197265625)};
    // // 6-point upwind biased
    // constexpr Array1D<Real,-2,3> d = {Real(-0.0093994140625), Real(0.0640869140652),
    //                                   Real(-0.201416015625), Real(0.469970703125),
    //                                   Real(0.704956054687), Real(-0.0281982421875)};
    // constexpr Array1D<Real,-2,3> c = {Real(-0.0076904296875), Real(0.0555419921875),
    //                                   Real(-0.199951171875), Real(0.999755859375),
    //                                   Real(0.1666259765625), Real(-0.0142822265625)}; 
    // // 7-point upwind biased
    // constexpr Array1D<Real,-3,3> d = {Real(-0.003524780273438), Real(0.028839111328125),
    //                                   Real(-0.113296508789063), Real(0.352478027343750),
    //                                   Real(0.793075561523438), Real(-0.063446044921875),
    //                                   Real(0.005874633789063)};
    // constexpr Array1D<Real,-3,3> c = {Real(-0.003524780273438), Real(0.030548095703125),
    //                                   Real(-0.137466430664063), Real(0.916442871093750),
    //                                   Real(0.229110717773438), Real(-0.039276123046875),
    //                                   Real(0.004165649414063)};
    
    // 6-point centered
    constexpr Array1D<Real,-3,2> c = {Real(0.007690429687500), Real(-0.06042480468750),
                                      Real(0.281982421875000), Real(0.845947265625000),
                                      Real(-0.08459472656250), Real(0.009399414062500)};    
    
    // 检查系数和是否为1
    // Real sum_c = 0.0;
    // for (int i = -3; i <= 2; ++i) {
    //     sum_c += c(i);
    // }
    // AMREX_ASSERT(std::abs(sum_c - 1.0) < 1e-10);
    // Real sum_d = 0.0;
    // for (int i = -3; i <= 2; ++i) {
    //     sum_d += d(i);
    // }
    // AMREX_ASSERT(std::abs(sum_d - 1.0) < 1e-10);
    // end checking

    int ii = amrex::coarsen(i,2);
    int s = i-ii*2-1;  // if i == ii*2, s = -1; if i == ii*2+1, s = 0;
    if (s == -1)
        // // overall 5-point upwind, 1-point biased
        // fine(i,j,k,n) = d(-2)*crse(ii-3,j,k,n)
        //     +           d(-1)*crse(ii-2,j,k,n)
        //     +           d( 0)*crse(ii-1,j,k,n)
        //     +           d( 1)*crse(ii  ,j,k,n)
        //     +           d( 2)*crse(ii+1,j,k,n);

        // // overall 6-point upwind, 2-point biased
        // fine(i,j,k,n) = d(-2)*crse(ii-4,j,k,n)
        //     +           d(-1)*crse(ii-3,j,k,n)
        //     +           d( 0)*crse(ii-2,j,k,n)
        //     +           d( 1)*crse(ii-1,j,k,n)
        //     +           d( 2)*crse(ii  ,j,k,n)
        //     +           d( 3)*crse(ii+1,j,k,n);

        // // overall 7-point upwind, 1-point biased
        // fine(i,j,k,n) = d(-3)*crse(ii-4,j,k,n)
        //     +           d(-2)*crse(ii-3,j,k,n)
        //     +           d(-1)*crse(ii-2,j,k,n)
        //     +           d( 0)*crse(ii-1,j,k,n)
        //     +           d( 1)*crse(ii  ,j,k,n)
        //     +           d( 2)*crse(ii+1,j,k,n)
        //     +           d( 3)*crse(ii+2,j,k,n);   

        // 6-point centered
        fine(i,j,k,n) = c(-3)*crse(ii-3,j,k,n)
            +           c(-2)*crse(ii-2,j,k,n)
            +           c(-1)*crse(ii-1,j,k,n)
            +           c( 0)*crse(ii  ,j,k,n)
            +           c( 1)*crse(ii+1,j,k,n)
            +           c( 2)*crse(ii+2,j,k,n);
    else if (s == 0)
        // // overall 5-point upwind, 1-point biased
        // fine(i,j,k,n) = c(-2)*crse(ii-2,j,k,n)
        //     +           c(-1)*crse(ii-1,j,k,n)
        //     +           c( 0)*crse(ii  ,j,k,n)
        //     +           c( 1)*crse(ii+1,j,k,n)
        //     +           c( 2)*crse(ii+2,j,k,n);

        // // overall 6-point upwind, 2-point biased
        // fine(i,j,k,n) = c(-2)*crse(ii-3,j,k,n)
        //     +           c(-1)*crse(ii-2,j,k,n)
        //     +           c( 0)*crse(ii-1,j,k,n)
        //     +           c( 1)*crse(ii  ,j,k,n)
        //     +           c( 2)*crse(ii+1,j,k,n)
        //     +           c( 3)*crse(ii+2,j,k,n);

        // // overall 7-point upwind, 1-point biased
        // fine(i,j,k,n) = c(-3)*crse(ii-3,j,k,n)
        //     +           c(-2)*crse(ii-2,j,k,n)
        //     +           c(-1)*crse(ii-1,j,k,n)
        //     +           c( 0)*crse(ii  ,j,k,n)
        //     +           c( 1)*crse(ii+1,j,k,n)
        //     +           c( 2)*crse(ii+2,j,k,n)
        //     +           c( 3)*crse(ii+3,j,k,n); 

        //  6-point centered
        fine(i,j,k,n) = c(-3)*crse(ii+3,j,k,n)
            +           c(-2)*crse(ii+2,j,k,n)
            +           c(-1)*crse(ii+1,j,k,n)
            +           c( 0)*crse(ii  ,j,k,n)
            +           c( 1)*crse(ii-1,j,k,n)
            +           c( 2)*crse(ii-2,j,k,n);
    else
        amrex::Abort("Lagrange5_interp: how did this happen?");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Lagrange5_interp_y (int i, int j, int k, int n, Array4<Real> const& fine,
                         Array4<Real const> const& crse) noexcept
{
    // 6-point centered
    constexpr Array1D<Real,-3,2> c = {Real(0.007690429687500), Real(-0.06042480468750),
                                      Real(0.281982421875000), Real(0.845947265625000),
                                      Real(-0.08459472656250), Real(0.009399414062500)}; 
    int jj = amrex::coarsen(j,2);
    int s = j-jj*2-1;  // if j == jj*2, s = -1; if j == jj*2+1, s = 0;
    if (s == -1)
        fine(i,j,k,n) = c(-3)*crse(i,jj-3,k,n)
            +           c(-2)*crse(i,jj-2,k,n)
            +           c(-1)*crse(i,jj-1,k,n)
            +           c( 0)*crse(i,jj  ,k,n)
            +           c( 1)*crse(i,jj+1,k,n)
            +           c( 2)*crse(i,jj+2,k,n);
    else if (s == 0)
        fine(i,j,k,n) = c(-3)*crse(i,jj+3,k,n)
            +           c(-2)*crse(i,jj+2,k,n)
            +           c(-1)*crse(i,jj+1,k,n)
            +           c( 0)*crse(i,jj  ,k,n)
            +           c( 1)*crse(i,jj+1,k,n)
            +           c( 2)*crse(i,jj+2,k,n);
    else
        amrex::Abort("Lagrange5_interp: how did this happen?");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Lagrange5_interp_z (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    // 6-point centered
    constexpr Array1D<Real,-3,2> c = {Real(0.007690429687500), Real(-0.06042480468750),
                                      Real(0.281982421875000), Real(0.845947265625000),
                                      Real(-0.08459472656250), Real(0.009399414062500)}; 
    int kk = amrex::coarsen(k,2);
    int s = k-kk*2-1;  // if k == kk*2, s = -1; if k == kk*2+1, s = 0;
    if (s == 0)
        fine(i,j,k,n) = c(-3)*crse(i,j,kk-3,n)
            +           c(-2)*crse(i,j,kk-2,n)
            +           c(-1)*crse(i,j,kk-1,n)
            +           c( 0)*crse(i,j,k   ,n)
            +           c( 1)*crse(i,j,kk+1,n)
            +           c( 2)*crse(i,j,kk+2,n);
    else if (s == -1)
        fine(i,j,k,n) = c(-3)*crse(i,j,kk+3,n)
            +           c(-2)*crse(i,j,kk+2,n)
            +           c(-1)*crse(i,j,kk+1,n)
            +           c( 0)*crse(i,j,k   ,n)
            +           c( 1)*crse(i,j,kk-1,n)
            +           c( 2)*crse(i,j,kk-2,n);
    else
        amrex::Abort("Lagrange5_interp: how did this happen?");
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_restrict_x (int i, int j, int k, int n, Array4<Real> const& crse,
                            Array4<Real const> const& fine) noexcept
{
    constexpr Array1D<Real,-2,3> c = {Real(0.01171875),  Real(-0.09765625),
                                      Real(0.5859375),   Real(0.5859375),
                                      Real(-0.09765625), Real(0.01171875)};
    constexpr Array1D<Real,-2,2> b1 = {Real(0.0234375), Real(-0.15625),
                                      Real(0.703125), Real(0.46875),
                                      Real(-0.0390625)}; // 5-point stencil
    constexpr Array1D<Real,-2,3> b2 = {Real(0.24609375), Real(1.23046875),
                                        Real(-0.8203125), Real(0.4921875),
                                        Real(-0.17578125), Real(0.02734375)};   // 6-point stencil    
    constexpr Array1D<Real,-3,2> d = {Real(-0.0093994140625), Real(0.0640869140652),
                                      Real(-0.201416015625), Real(0.469970703125),
                                      Real(0.704956054687), Real(-0.0281982421875)}; // 6-point, biased stencil
    constexpr Array1D<Real,-3,4> c8 = {Real(-0.002441406250000), Real(0.023925781250000),
                                      Real(-0.119628906250000), Real(0.598144531250000),
                                      Real(0.598144531250000), Real(-0.119628906250000),
                                      Real(0.023925781250000), Real(-0.002441406250000)}; // 8-point, centered stencil

    int small_end = crse.begin.x;
    int big_end   = crse.end.x - 1;
    if (i > small_end && i < big_end) {
        int ii = i*2;
        // 8-point centered
        // crse(i,j,k,n) = c8(-3)*fine(ii-3,j,k,n)
        //     +           c8(-2)*fine(ii-2,j,k,n)
        //     +           c8(-1)*fine(ii-1,j,k,n)
        //     +           c8( 0)*fine(ii  ,j,k,n)
        //     +           c8( 1)*fine(ii+1,j,k,n)
        //     +           c8( 2)*fine(ii+2,j,k,n)
        //     +           c8( 3)*fine(ii+3,j,k,n);
        //     +           c8( 4)*fine(ii+4,j,k,n);
        // 6-point centered
        crse(i,j,k,n) = c(-2)*fine(ii-2,j,k,n)
            +           c(-1)*fine(ii-1,j,k,n)
            +           c( 0)*fine(ii  ,j,k,n)
            +           c( 1)*fine(ii+1,j,k,n)
            +           c( 2)*fine(ii+2,j,k,n)
            +           c( 3)*fine(ii+3,j,k,n);    
        // // 6-point biased
        // crse(i,j,k,n) = d(-3)*fine(ii-3,j,k,n)
        //     +           d(-2)*fine(ii-2,j,k,n)
        //     +           d(-1)*fine(ii-1,j,k,n)
        //     +           d( 0)*fine(ii  ,j,k,n)
        //     +           d( 1)*fine(ii+1,j,k,n)
        //     +           d( 2)*fine(ii+2,j,k,n);    
        // 5-point biased
        // crse(i,j,k,n) = b1(-2)*fine(ii-2,j,k,n)
        //     +           b1(-1)*fine(ii-1,j,k,n)
        //     +           b1( 0)*fine(ii  ,j,k,n)
        //     +           b1( 1)*fine(ii+1,j,k,n)
        //     +           b1( 2)*fine(ii+2,j,k,n);    
    } else if (i == small_end) {
        int ii = i*2 + 2;
        crse(i,j,k,n) = b2(-2)*fine(ii-2,j,k,n)
            +           b2(-1)*fine(ii-1,j,k,n)
            +           b2( 0)*fine(ii  ,j,k,n)
            +           b2( 1)*fine(ii+1,j,k,n)
            +           b2( 2)*fine(ii+2,j,k,n)
            +           b2( 3)*fine(ii+3,j,k,n);
    } else if (i == big_end) {
        int ii = i*2 - 1;
        crse(i,j,k,n) = b2(-2)*fine(ii+2,j,k,n)
            +           b2(-1)*fine(ii+1,j,k,n)
            +           b2( 0)*fine(ii  ,j,k,n)
            +           b2( 1)*fine(ii-1,j,k,n)
            +           b2( 2)*fine(ii-2,j,k,n)
            +           b2( 3)*fine(ii-3,j,k,n);
    } else {
        amrex::Abort("cell_quartic_restrict_x:  how did this happen?");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_restrict_y (int i, int j, int k, int n, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    constexpr Array1D<Real,-2,3> c  = {Real(0.01171875),  Real(-0.09765625),
                                       Real(0.5859375),   Real(0.5859375),
                                       Real(-0.09765625), Real(0.01171875)};
    constexpr Array1D<Real,-2,2> b1 = {Real(0.2734375), Real(1.09375),
                                       Real(-0.546875), Real(0.21875),
                                       Real(-0.0390625)}; // 5-point stencil
    constexpr Array1D<Real,-2,3> b2 = {Real(0.24609375), Real(1.23046875),
                                       Real(-0.8203125), Real(0.4921875),
                                       Real(-0.17578125), Real(0.02734375)};   // 6-point stencil    

    int small_end = crse.begin.y, big_end = crse.end.y - 1;
    if (j > small_end && j < big_end) {
        int jj = j*2;
        crse(i,j,k,n) = c(-2)*fine(i,jj-2,k,n)
            +           c(-1)*fine(i,jj-1,k,n)
            +           c( 0)*fine(i,jj,  k,n)
            +           c( 1)*fine(i,jj+1,k,n)
            +           c( 2)*fine(i,jj+2,k,n)
            +           c( 3)*fine(i,jj+3,k,n);
    } else if (j == small_end) {
        int jj = j*2 + 2;
        crse(i,j,k,n) = b2(-2)*fine(i,jj-2,k,n)
            +           b2(-1)*fine(i,jj-1,k,n)
            +           b2( 0)*fine(i,jj,  k,n)
            +           b2( 1)*fine(i,jj+1,k,n)
            +           b2( 2)*fine(i,jj+2,k,n)
            +           b2( 3)*fine(i,jj+3,k,n);
    } else if (j == big_end) {
        int jj = j*2 - 1;
        crse(i,j,k,n) = b2(-2)*fine(i,jj+2,k,n)
            +           b2(-1)*fine(i,jj+1,k,n)
            +           b2( 0)*fine(i,jj,  k,n)
            +           b2( 1)*fine(i,jj-1,k,n)
            +           b2( 2)*fine(i,jj-2,k,n)
            +           b2( 3)*fine(i,jj-3,k,n);
    } else {
        amrex::Abort("cell_quartic_restrict_y:  how did this happen?");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_restrict_z (int i, int j, int k, int n, Array4<Real> const& crse,
                            Array4<Real const> const& fine) noexcept
{
    constexpr Array1D<Real,-2,3> c = {Real(0.01171875),  Real(-0.09765625),
                                      Real(0.5859375),   Real(0.5859375),
                                      Real(-0.09765625), Real(0.01171875)};
    constexpr Array1D<Real,-2,2> b1 = {Real(0.2734375), Real(1.09375),
                                      Real(-0.546875), Real(0.21875),
                                      Real(-0.0390625)}; // 5-point stencil
    constexpr Array1D<Real,-2,3> b2 = {Real(0.24609375), Real(1.23046875),
                                        Real(-0.8203125), Real(0.4921875),
                                        Real(-0.17578125), Real(0.02734375)};   // 6-point stencil    

    int small_end = crse.begin.z, big_end = crse.end.z - 1;
    if (k > small_end && k < big_end) {
        int kk = k*2;
        crse(i,j,k,n) = c(-2)*fine(i,j,kk-2,n)
            +           c(-1)*fine(i,j,kk-1,n)
            +           c( 0)*fine(i,j,kk  ,n)
            +           c( 1)*fine(i,j,kk+1,n)
            +           c( 2)*fine(i,j,kk+2,n)
            +           c( 3)*fine(i,j,kk+3,n);
    } else if (k == small_end) {
        int kk = k*2 + 2;
        crse(i,j,k,n) = b2(-2)*fine(i,j,kk-2,n)
            +           b2(-1)*fine(i,j,kk-1,n)
            +           b2( 0)*fine(i,j,kk  ,n)
            +           b2( 1)*fine(i,j,kk+1,n)
            +           b2( 2)*fine(i,j,kk+2,n)
            +           b2( 3)*fine(i,j,kk+3,n);
    } else if (k == big_end) {
        int kk = k*2 - 1;
        crse(i,j,k,n) = b2(-2)*fine(i,j,kk+2,n)
            +           b2(-1)*fine(i,j,kk+1,n)
            +           b2( 0)*fine(i,j,kk  ,n)
            +           b2( 1)*fine(i,j,kk-1,n)
            +           b2( 2)*fine(i,j,kk-2,n)
            +           b2( 3)*fine(i,j,kk-3,n);
    } else {
        amrex::Abort("cell_quartic_restrict_z:  how did this happen?");
    }
}

// WENO coefficients stored in constant memory for optimal GPU performance
namespace {
#ifdef AMREX_USE_GPU
    // Linear weights for different ratio and sub-cell positions
    // [ratio_idx][subcell_idx][stencil_idx]
    AMREX_GPU_CONSTANT Real weno_linear_weights[2][4][4] = {
        // ratio = 2
        {
            {3.0/256.0, 99.0/256.0, 693.0/1280.0, 77.0/1280.0},  // subcell 0 (r=2, s=-1)
            {3.0/256.0, 99.0/256.0, 693.0/1280.0, 77.0/1280.0},  // subcell 1 (r=2, s=1)
            {0.0, 0.0, 0.0, 0.0},  // unused
            {0.0, 0.0, 0.0, 0.0}   // unused
        },
        // ratio = 4
        {
            {209.0/10240.0, 1462.0/3411.0, 271.0/535.0, 91.0/2480.0},  // subcell 0 (r=4, ii)
            {5.0/1004.0, 111.0/323.0, 1173.0/2048.0, 68.0/865.0},      // subcell 1 (r=4, ii+1/2)
            {5.0/1004.0, 111.0/323.0, 1173.0/2048.0, 68.0/865.0},      // subcell 2 (r=4, ii+1)
            {209.0/10240.0, 1462.0/3411.0, 271.0/535.0, 91.0/2480.0}   // subcell 3 (r=4, ii+3/2)
        }
    };

    // Smoothness indicator coefficients [coeff_idx][term_idx]
    // beta = c0*v0*v2 + c1*v0*v1 + c2*v1*v2 + c3*v0^2 + c4*v1^2 + c5*v2^2
    AMREX_GPU_CONSTANT Real weno_beta_coeffs[4][6] = {
        {29.0/3.0, -49.0/3.0, -73.0/3.0, 10.0/3.0, 61.0/3.0, 22.0/3.0},  // beta0
        {11.0/3.0, -19.0/3.0, -31.0/3.0,  4.0/3.0, 25.0/3.0, 10.0/3.0},  // beta1
        { 5.0/3.0, -13.0/3.0, -13.0/3.0,  4.0/3.0, 13.0/3.0,  4.0/3.0},  // beta2
        {11.0/3.0, -31.0/3.0, -19.0/3.0, 10.0/3.0, 25.0/3.0,  4.0/3.0}   // beta3
    };

    // Stencil coefficients for different ratio and sub-cell positions
    // [ratio_idx][subcell_idx][stencil_idx][coeff_idx]
    AMREX_GPU_CONSTANT Real weno_stencil_coeffs[2][4][4][3] = {
        // ratio = 2
        {
            // subcell 0 (r=2, s=-1)
            {{21.0/32.0, -66.0/32.0, 77.0/32.0}, {-3.0/32.0, 14.0/32.0, 21.0/32.0}, 
             {5.0/32.0, 30.0/32.0, -3.0/32.0}, {45.0/32.0, -18.0/32.0, 5.0/32.0}},
            // subcell 1 (r=2, s=1) - same as subcell 0 for symmetric case
            {{21.0/32.0, -66.0/32.0, 77.0/32.0}, {-3.0/32.0, 14.0/32.0, 21.0/32.0}, 
             {5.0/32.0, 30.0/32.0, -3.0/32.0}, {45.0/32.0, -18.0/32.0, 5.0/32.0}},
            // unused
            {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
            {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
        },
        // ratio = 4
        {
            // subcell 0 (r=4, ii)
            {{65.0/128.0, -210.0/128.0, 273.0/128.0}, {-15.0/128.0, 78.0/128.0, 65.0/128.0}, 
             {33.0/128.0, 110.0/128.0, -15.0/128.0}, {209.0/128.0, -114.0/128.0, 33.0/128.0}},
            // subcell 1 (r=4, ii+1/2)
            {{105.0/128.0, -322.0/128.0, 345.0/128.0}, {-7.0/128.0, 30.0/128.0, 105.0/128.0}, 
             {9.0/128.0, 126.0/128.0, -7.0/128.0}, {153.0/128.0, -34.0/128.0, 9.0/128.0}},
            // subcell 2 (r=4, ii+1) - mirror of subcell 1
            {{105.0/128.0, -322.0/128.0, 345.0/128.0}, {-7.0/128.0, 30.0/128.0, 105.0/128.0}, 
             {9.0/128.0, 126.0/128.0, -7.0/128.0}, {153.0/128.0, -34.0/128.0, 9.0/128.0}},
            // subcell 3 (r=4, ii+3/2) - mirror of subcell 0
            {{65.0/128.0, -210.0/128.0, 273.0/128.0}, {-15.0/128.0, 78.0/128.0, 65.0/128.0}, 
             {33.0/128.0, 110.0/128.0, -15.0/128.0}, {209.0/128.0, -114.0/128.0, 33.0/128.0}}
        }
    };
#else
    // Linear weights for different ratio and sub-cell positions
    // [ratio_idx][subcell_idx][stencil_idx]
    Real weno_linear_weights[2][4][4] = {
        // ratio = 2
        {
            {3.0/256.0, 99.0/256.0, 693.0/1280.0, 77.0/1280.0},  // subcell 0 (r=2, s=-1)
            {3.0/256.0, 99.0/256.0, 693.0/1280.0, 77.0/1280.0},  // subcell 1 (r=2, s=1)
            {0.0, 0.0, 0.0, 0.0},  // unused
            {0.0, 0.0, 0.0, 0.0}   // unused
        },
        // ratio = 4
        {
            {209.0/10240.0, 1462.0/3411.0, 271.0/535.0, 91.0/2480.0},  // subcell 0 (r=4, ii)
            {5.0/1004.0, 111.0/323.0, 1173.0/2048.0, 68.0/865.0},      // subcell 1 (r=4, ii+1/2)
            {5.0/1004.0, 111.0/323.0, 1173.0/2048.0, 68.0/865.0},      // subcell 2 (r=4, ii+1)
            {209.0/10240.0, 1462.0/3411.0, 271.0/535.0, 91.0/2480.0}   // subcell 3 (r=4, ii+3/2)
        }
    };

    // Smoothness indicator coefficients [coeff_idx][term_idx]
    // beta = c0*v0*v2 + c1*v0*v1 + c2*v1*v2 + c3*v0^2 + c4*v1^2 + c5*v2^2
    Real weno_beta_coeffs[4][6] = {
        {29.0/3.0, -49.0/3.0, -73.0/3.0, 10.0/3.0, 61.0/3.0, 22.0/3.0},  // beta0
        {11.0/3.0, -19.0/3.0, -31.0/3.0,  4.0/3.0, 25.0/3.0, 10.0/3.0},  // beta1
        { 5.0/3.0, -13.0/3.0, -13.0/3.0,  4.0/3.0, 13.0/3.0,  4.0/3.0},  // beta2
        {11.0/3.0, -31.0/3.0, -19.0/3.0, 10.0/3.0, 25.0/3.0,  4.0/3.0}   // beta3
    };

    // Stencil coefficients for different ratio and sub-cell positions
    // [ratio_idx][subcell_idx][stencil_idx][coeff_idx]
    Real weno_stencil_coeffs[2][4][4][3] = {
        // ratio = 2
        {
            // subcell 0 (r=2, s=-1)
            {{21.0/32.0, -66.0/32.0, 77.0/32.0}, {-3.0/32.0, 14.0/32.0, 21.0/32.0}, 
             {5.0/32.0, 30.0/32.0, -3.0/32.0}, {45.0/32.0, -18.0/32.0, 5.0/32.0}},
            // subcell 1 (r=2, s=1) - same as subcell 0 for symmetric case
            {{21.0/32.0, -66.0/32.0, 77.0/32.0}, {-3.0/32.0, 14.0/32.0, 21.0/32.0}, 
             {5.0/32.0, 30.0/32.0, -3.0/32.0}, {45.0/32.0, -18.0/32.0, 5.0/32.0}},
            // unused
            {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
            {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
        },
        // ratio = 4
        {
            // subcell 0 (r=4, ii)
            {{65.0/128.0, -210.0/128.0, 273.0/128.0}, {-15.0/128.0, 78.0/128.0, 65.0/128.0}, 
             {33.0/128.0, 110.0/128.0, -15.0/128.0}, {209.0/128.0, -114.0/128.0, 33.0/128.0}},
            // subcell 1 (r=4, ii+1/2)
            {{105.0/128.0, -322.0/128.0, 345.0/128.0}, {-7.0/128.0, 30.0/128.0, 105.0/128.0}, 
             {9.0/128.0, 126.0/128.0, -7.0/128.0}, {153.0/128.0, -34.0/128.0, 9.0/128.0}},
            // subcell 2 (r=4, ii+1) - mirror of subcell 1
            {{105.0/128.0, -322.0/128.0, 345.0/128.0}, {-7.0/128.0, 30.0/128.0, 105.0/128.0}, 
             {9.0/128.0, 126.0/128.0, -7.0/128.0}, {153.0/128.0, -34.0/128.0, 9.0/128.0}},
            // subcell 3 (r=4, ii+3/2) - mirror of subcell 0
            {{65.0/128.0, -210.0/128.0, 273.0/128.0}, {-15.0/128.0, 78.0/128.0, 65.0/128.0}, 
             {33.0/128.0, 110.0/128.0, -15.0/128.0}, {209.0/128.0, -114.0/128.0, 33.0/128.0}}
        }
    };
#endif
}

// Unified WENO interpolation function supporting both ratio=2 and ratio=4
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void WENO_Interp_Unified(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
                         const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
                         int ratio_idx, int subcell_idx, amrex::Real& v_edge)
{
    using namespace amrex;
    
    // Compute smoothness indicators using constant memory coefficients
    Real beta[4];
    Real tvals[6] = {v0, v1, v2, v3, v4, v5};
    
    for (int k = 0; k < 4; ++k) {
        const Real* coeffs = weno_beta_coeffs[k];
        Real vals[3] = {tvals[k], tvals[k+1], tvals[k+2]}; // Center values for beta calculation
        beta[k] = coeffs[0]*vals[0]*vals[2] + coeffs[1]*vals[0]*vals[1] + coeffs[2]*vals[1]*vals[2] +
                  coeffs[3]*vals[0]*vals[0] + coeffs[4]*vals[1]*vals[1] + coeffs[5]*vals[2]*vals[2];
    }
    
    // Small value for computing non-linear weights
    constexpr Real eps = 1.0e-6;
    
    // Compute alpha values and normalize
    Real alpha[4], sum_alpha = 0.0;
    for (int k = 0; k < 4; ++k) {
        alpha[k] = weno_linear_weights[ratio_idx][subcell_idx][k] / ((eps + beta[k]) * (eps + beta[k]));
        sum_alpha += alpha[k];
    }
    
    // Compute candidate stencils and final result
    Real result = 0.0;
    
    for (int k = 0; k < 4; ++k) {
        Real omega = alpha[k] / sum_alpha;
        Real q = 0.0;
        for (int j = 0; j < 3; ++j) {
            q += weno_stencil_coeffs[ratio_idx][subcell_idx][k][j] * tvals[k + j];
        }
        result += omega * q;
    }
    
    v_edge = result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_x (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse,
                            IntVect const& ratio) noexcept
{
    if (ratio[0] == 2) {
        // ratio = 2 case
        int ii = amrex::coarsen(i,2);
        int s = 2 * (i - ii*2) - 1; // i=ii*2, s = -1; i=ii*2+1, s = 1
        int subcell_idx = i - ii*2;

        WENO_Interp_Unified(crse(ii+s*3,j,k,n), crse(ii+s*2,j,k,n), 
                           crse(ii+s*1,j,k,n), crse(ii,j,k,n),
                           crse(ii-s*1,j,k,n), crse(ii-s*2,j,k,n),
                           0, subcell_idx, fine(i,j,k,n));
    } else if (ratio[0] == 4) {
        // ratio = 4 case
        int ii = amrex::coarsen(i,4);
        int subcell_idx = i - ii*4; // 0, 1, 2, or 3
        int s = 2 * (subcell_idx/2) - 1; // i=ii*2 or ii*2+1, s = -1; i=ii*2+2 or ii*2+3, s = 1
        
        WENO_Interp_Unified(crse(ii+s*3,j,k,n), crse(ii+s*2,j,k,n), 
                           crse(ii+s*1,j,k,n), crse(ii,j,k,n),
                           crse(ii-s*1,j,k,n), crse(ii-s*2,j,k,n),
                           1, subcell_idx, fine(i,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_y (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse,
                            IntVect const& ratio) noexcept
{
    if (ratio[1] == 2) {
        // ratio = 2 case
        int jj = amrex::coarsen(j,2);
        int s = 2 * (j - jj*2) - 1; // j=jj*2, s = -1; j=jj*2+1, s = 1
        int subcell_idx = j - jj*2;

        WENO_Interp_Unified(crse(i,jj+s*3,k,n), crse(i,jj+s*2,k,n), 
                           crse(i,jj+s*1,k,n), crse(i,jj,k,n),
                           crse(i,jj-s*1,k,n), crse(i,jj-s*2,k,n),
                           0, subcell_idx, fine(i,j,k,n));
    } else if (ratio[1] == 4) {
        // ratio = 4 case
        int jj = amrex::coarsen(j,4);
        int subcell_idx = j - jj*4; // 0, 1, 2, or 3
        int s = 2 * (subcell_idx/2) - 1; // j=jj*2 or jj*2+1, s = -1; j=jj*2+2 or jj*2+3, s = 1
        
        WENO_Interp_Unified(crse(i,jj+s*3,k,n), crse(i,jj+s*2,k,n), 
                           crse(i,jj+s*1,k,n), crse(i,jj,k,n),
                           crse(i,jj-s*1,k,n), crse(i,jj-s*2,k,n),
                           1, subcell_idx, fine(i,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_z (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse,
                            IntVect const& ratio) noexcept
{
    if (ratio[2] == 2) {
        // ratio = 2 case
        int kk = amrex::coarsen(k,2);
        int s = 2 * (k - kk*2) - 1; // k=kk*2, s = -1; k=kk*2+1, s = 1
        int subcell_idx = k - kk*2;

        WENO_Interp_Unified(crse(i,j,kk+s*3,n), crse(i,j,kk+s*2,n), 
                           crse(i,j,kk+s*1,n), crse(i,j,kk,n),
                           crse(i,j,kk-s*1,n), crse(i,j,kk-s*2,n),
                           0, subcell_idx, fine(i,j,k,n));
    } else if (ratio[2] == 4) {
        // ratio = 4 case
        int kk = amrex::coarsen(k,4);
        int subcell_idx = k - kk*4; // 0, 1, 2, or 3
        int s = 2 * (subcell_idx/2) - 1; // k=kk*2 or kk*2+1, s = -1; k=kk*2+2 or kk*2+3, s = 1
        
        WENO_Interp_Unified(crse(i,j,kk+s*3,n), crse(i,j,kk+s*2,n), 
                           crse(i,j,kk+s*1,n), crse(i,j,kk,n),
                           crse(i,j,kk-s*1,n), crse(i,j,kk-s*2,n),
                           1, subcell_idx, fine(i,j,k,n));
    }
}

// WENO restriction coefficients stored in constant memory
namespace {
#ifdef AMREX_USE_GPU
    // Linear weights for restriction operations [ratio_idx][use_eno]
    AMREX_GPU_CONSTANT Real weno_restrict_linear_weights[2][2][4] = {
        // ratio = 2
        {
            // central WENO scheme
            {1.0/32.0, 15.0/32.0, 15.0/32.0, 1.0/32.0},
            // boundary WENO scheme  
            {21.0/32.0, 63.0/160.0, -9.0/160.0, 1.0/160.0}
        },
        // ratio = 4
        {
            // central WENO scheme
            {1.0/32.0, 15.0/32.0, 15.0/32.0, 1.0/32.0},
            // boundary WENO scheme
            {7.0/32.0, 21.0/32.0, 21.0/160.0, -1.0/160.0}
        }
    };

    // Smoothness indicator coefficients [coeff_idx][term_idx]
    // beta = c0*v0*v2 + c1*v0*v1 + c2*v1*v2 + c3*v0^2 + c4*v1^2 + c5*v2^2
    AMREX_GPU_CONSTANT Real weno_restrict_beta_coeffs[2][4][6] = {
        // central
        {
            {74.0/12.0, -124.0/12.0, -196.0/12.0, 25.0/12.0, 160.0/12.0, 61.0/12.0},  // beta0
            {26.0/12.0,  -52.0/12.0,  -76.0/12.0, 13.0/12.0,  64.0/12.0, 25.0/12.0},  // beta1
            {26.0/12.0,  -76.0/12.0,  -52.0/12.0, 25.0/12.0,  64.0/12.0, 13.0/12.0},  // beta2
            {74.0/12.0, -196.0/12.0, -124.0/12.0, 61.0/12.0, 160.0/12.0, 25.0/12.0}   // beta3
        },
        // boundary
        {
            {26.0/12.0,   -52.0/12.0,  -76.0/12.0,  13.0/12.0,  64.0/12.0, 25.0/12.0},
            {26.0/12.0,   -76.0/12.0,  -52.0/12.0,  25.0/12.0,  64.0/12.0, 13.0/12.0},
            {74.0/12.0,  -196.0/12.0, -124.0/12.0,  61.0/12.0, 160.0/12.0, 25.0/12.0},
            {170.0/12.0, -412.0/12.0, -292.0/12.0, 121.0/12.0, 352.0/12.0, 61.0/12.0}
        }
    };

    // Stencil coefficients for restriction [ratio_idx][use_eno][stencil_idx][coeff_idx]
    AMREX_GPU_CONSTANT Real weno_restrict_stencil_coeffs[2][2][4][3] = {
        // ratio = 2
        {
            // central WENO scheme
            {
                {3.0/8.0, -10.0/8.0, 15.0/8.0},    // q0
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q1  
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q2
                {15.0/8.0, -10.0/8.0, 3.0/8.0}     // q3
            },
            // boundary WENO scheme
            {
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q0
                {15.0/8.0, -10.0/8.0, 3.0/8.0},    // q1  
                {35.0/8.0, -42.0/8.0, 15.0/8.0},   // q2
                {63.0/8.0, -90.0/8.0, 35.0/8.0}    // q3
            }
        },
        // ratio = 4
        {
            // central WENO scheme
            {
                {3.0/8.0, -10.0/8.0, 15.0/8.0},    // q0
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q1
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q2
                {15.0/8.0, -10.0/8.0, 3.0/8.0}     // q3
            },
            // boundary WENO scheme
            {
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q0
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q1
                {15.0/8.0, -10.0/8.0, 3.0/8.0},    // q2
                {35.0/8.0, -42.0/8.0, 15.0/8.0}    // q3
            }
        }
    };

    // Simple 2nd order Lagrange coefficients for boundary cases [ratio_idx]
    AMREX_GPU_CONSTANT Real eno_restrict_coeffs[2][4] = {
        // // ratio = 2
        // {3.0/8.0, 6.0/8.0, -1.0/8.0, 0.0},
        // // ratio = 4
        // {-1.0/8.0, 6.0/8.0, 3.0/8.0, 0.0}
        // ratio = 2
        {1.0/2.0, 1.0/2.0, 0.0, 0.0},
        // ratio = 4
        {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0}
    };
#else
    // Linear weights for restriction operations [ratio_idx][use_eno]
    Real weno_restrict_linear_weights[2][2][4] = {
        // ratio = 2
        {
            // central WENO scheme
            {1.0/32.0, 15.0/32.0, 15.0/32.0, 1.0/32.0},
            // boundary WENO scheme  
            {21.0/32.0, 63.0/160.0, -9.0/160.0, 1.0/160.0}
        },
        // ratio = 4
        {
            // central WENO scheme
            {1.0/32.0, 15.0/32.0, 15.0/32.0, 1.0/32.0},
            // boundary WENO scheme
            {7.0/32.0, 21.0/32.0, 21.0/160.0, -1.0/160.0}
        }
    };

    // Smoothness indicator coefficients [coeff_idx][term_idx]
    // beta = c0*v0*v2 + c1*v0*v1 + c2*v1*v2 + c3*v0^2 + c4*v1^2 + c5*v2^2
    Real weno_restrict_beta_coeffs[2][4][6] = {
        // central
        {
            {74.0/12.0, -124.0/12.0, -196.0/12.0, 25.0/12.0, 160.0/12.0, 61.0/12.0},  // beta0
            {26.0/12.0,  -52.0/12.0,  -76.0/12.0, 13.0/12.0,  64.0/12.0, 25.0/12.0},  // beta1
            {26.0/12.0,  -76.0/12.0,  -52.0/12.0, 25.0/12.0,  64.0/12.0, 13.0/12.0},  // beta2
            {74.0/12.0, -196.0/12.0, -124.0/12.0, 61.0/12.0, 160.0/12.0, 25.0/12.0}   // beta3
        },
        // boundary
        {
            {26.0/12.0,   -52.0/12.0,  -76.0/12.0,  13.0/12.0,  64.0/12.0, 25.0/12.0},
            {26.0/12.0,   -76.0/12.0,  -52.0/12.0,  25.0/12.0,  64.0/12.0, 13.0/12.0},
            {74.0/12.0,  -196.0/12.0, -124.0/12.0,  61.0/12.0, 160.0/12.0, 25.0/12.0},
            {170.0/12.0, -412.0/12.0, -292.0/12.0, 121.0/12.0, 352.0/12.0, 61.0/12.0}
        }
    };

    // Stencil coefficients for restriction [ratio_idx][use_eno][stencil_idx][coeff_idx]
    Real weno_restrict_stencil_coeffs[2][2][4][3] = {
        // ratio = 2
        {
            // central WENO scheme
            {
                {3.0/8.0, -10.0/8.0, 15.0/8.0},    // q0
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q1  
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q2
                {15.0/8.0, -10.0/8.0, 3.0/8.0}     // q3
            },
            // boundary WENO scheme
            {
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q0
                {15.0/8.0, -10.0/8.0, 3.0/8.0},    // q1  
                {35.0/8.0, -42.0/8.0, 15.0/8.0},   // q2
                {63.0/8.0, -90.0/8.0, 35.0/8.0}    // q3
            }
        },
        // ratio = 4
        {
            // central WENO scheme
            {
                {3.0/8.0, -10.0/8.0, 15.0/8.0},    // q0
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q1  
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q2
                {15.0/8.0, -10.0/8.0, 3.0/8.0}     // q3
            },
            // boundary WENO scheme
            {
                {-1.0/8.0, 6.0/8.0, 3.0/8.0},      // q0
                {3.0/8.0, 6.0/8.0, -1.0/8.0},      // q1
                {15.0/8.0, -10.0/8.0, 3.0/8.0},    // q2
                {35.0/8.0, -42.0/8.0, 15.0/8.0}    // q3
            }
        }
    };

    // Simple 2nd order Lagrange coefficients for boundary cases [ratio_idx]
    Real eno_restrict_coeffs[2][4] = {
        // // ratio = 2
        // {3.0/8.0, 6.0/8.0, -1.0/8.0, 0.0},
        // // ratio = 4
        // {-1.0/8.0, 6.0/8.0, 3.0/8.0, 0.0}
        // ratio = 2
        {1.0/2.0, 1.0/2.0, 0.0, 0.0},
        // ratio = 4
        {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0}
    };
#endif
}

// Unified WENO restriction function supporting both ratio=2 and ratio=4
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void WENO_Restrict_Unified(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
                           const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
                           int ratio_idx, bool use_eno, bool is_boundary, amrex::Real& v_edge)
{
    using namespace amrex;
    
    if (is_boundary && use_eno) {
        // Simple 2nd-order Lagrange interpolation for boundary
        const Real* coeffs = eno_restrict_coeffs[ratio_idx];
        v_edge = coeffs[0]*v0 + coeffs[1]*v1 + coeffs[2]*v2 + coeffs[3]*v3;
        return;
    }

    // WENO reconstruction
    Real vals[6] = {v0, v1, v2, v3, v4, v5};
        
    // Select appropriate linear weights based on use_eno
    int pos_idx = is_boundary ? 1 : 0;

    // Compute smoothness indicators using the same beta coefficients as in WENO_Interp_Unified
    Real beta[4];
    for (int k = 0; k < 4; ++k) {
        const Real* coeffs = weno_restrict_beta_coeffs[pos_idx][k];
        Real v[3] = {vals[k], vals[k+1], vals[k+2]};
        beta[k] = coeffs[0]*v[0]*v[2] + coeffs[1]*v[0]*v[1] + coeffs[2]*v[1]*v[2] +
                  coeffs[3]*v[0]*v[0] + coeffs[4]*v[1]*v[1] + coeffs[5]*v[2]*v[2];
    }
    
    // Small value for computing non-linear weights
    constexpr Real eps = 1.0e-6;

    // Compute alpha values and normalize
    Real alpha[4], sum_alpha = 0.0;
    for (int k = 0; k < 4; ++k) {
        alpha[k] = weno_restrict_linear_weights[ratio_idx][pos_idx][k] / ((eps + beta[k]) * (eps + beta[k]));
        sum_alpha += alpha[k];
    }
    
    // Compute candidate stencils and final result
    Real result = 0.0;
    for (int k = 0; k < 4; ++k) {
        Real omega = alpha[k] / sum_alpha;
        Real q = 0.0;
        for (int j = 0; j < 3; ++j) {
            q += weno_restrict_stencil_coeffs[ratio_idx][pos_idx][k][j] * vals[k + j];
        }
        result += omega * q;
    }
    
    v_edge = result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_x (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine, IntVect const& ratio) noexcept
{
    int small_end = crse.begin.x, big_end = crse.end.x-1;
    int ii = i*ratio[0];
    if (ratio[0] == 2) {
        if (i > small_end && i < big_end) {
            WENO_Restrict_Unified(fine(ii-2,j,k,n), fine(ii-1,j,k,n),
                                  fine(ii,j,k,n),   fine(ii+1,j,k,n),
                                  fine(ii+2,j,k,n), fine(ii+3,j,k,n),
                                  0, false, false, crse(i,j,k,n));
        } else if (i == small_end) {
            WENO_Restrict_Unified(fine(ii,j,k,n),   fine(ii+1,j,k,n),
                                  fine(ii+2,j,k,n), fine(ii+3,j,k,n),
                                  fine(ii+4,j,k,n), fine(ii+5,j,k,n),
                                  0, true, true, crse(i,j,k,n));
        } else if (i == big_end) {
            WENO_Restrict_Unified(fine(ii+1,j,k,n), fine(ii,j,k,n),
                                  fine(ii-1,j,k,n), fine(ii-2,j,k,n),
                                  fine(ii-3,j,k,n), fine(ii-4,j,k,n),
                                  0, true, true, crse(i,j,k,n));
        }
    } else if (ratio[0] == 4) {
        if (i > small_end && i < big_end) {
            WENO_Restrict_Unified(fine(ii-1,j,k,n), fine(ii,j,k,n),
                                  fine(ii+1,j,k,n), fine(ii+2,j,k,n),
                                  fine(ii+3,j,k,n), fine(ii+4,j,k,n),
                                  1, false, false, crse(i,j,k,n));
        } else if (i == small_end) {
            WENO_Restrict_Unified(fine(ii,j,k,n),   fine(ii+1,j,k,n),
                                  fine(ii+2,j,k,n), fine(ii+3,j,k,n),
                                  fine(ii+4,j,k,n), fine(ii+5,j,k,n),
                                  1, true, true, crse(i,j,k,n));
        } else if (i == big_end) {
            WENO_Restrict_Unified(fine(ii+3,j,k,n), fine(ii+2,j,k,n),
                                  fine(ii+1,j,k,n), fine(ii,j,k,n),
                                  fine(ii-1,j,k,n), fine(ii-2,j,k,n),
                                  1, true, true, crse(i,j,k,n));
        }
    } else {
        amrex::Abort("weno_restrict_x: how did this happen?");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_y (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine, IntVect const& ratio) noexcept
{
    int small_end = crse.begin.y, big_end = crse.end.y-1;
    int jj = j*ratio[1];
    if (ratio[1] == 2) {
        if (j > small_end && j < big_end) {
            WENO_Restrict_Unified(fine(i,jj-2,k,n), fine(i,jj-1,k,n),
                                  fine(i,jj,k,n),   fine(i,jj+1,k,n),
                                  fine(i,jj+2,k,n), fine(i,jj+3,k,n),
                                  0, false, false, crse(i,j,k,n));
        } else if (j == small_end) {
            WENO_Restrict_Unified(fine(i,jj,k,n),   fine(i,jj+1,k,n),
                                  fine(i,jj+2,k,n), fine(i,jj+3,k,n),
                                  fine(i,jj+4,k,n), fine(i,jj+5,k,n),
                                  0, true, true, crse(i,j,k,n));
        } else if (j == big_end) {
            WENO_Restrict_Unified(fine(i,jj+1,k,n), fine(i,jj,k,n),
                                  fine(i,jj-1,k,n), fine(i,jj-2,k,n),
                                  fine(i,jj-3,k,n), fine(i,jj-4,k,n),
                                  0, true, true, crse(i,j,k,n));
        }
    } else if (ratio[1] == 4) {
        if (j > small_end && j < big_end) {
            WENO_Restrict_Unified(fine(i,jj-1,k,n), fine(i,jj,k,n),
                                  fine(i,jj+1,k,n), fine(i,jj+2,k,n),
                                  fine(i,jj+3,k,n), fine(i,jj+4,k,n),
                                  1, false, false, crse(i,j,k,n));
        } else if (j == small_end) {
            WENO_Restrict_Unified(fine(i,jj,k,n),   fine(i,jj+1,k,n),
                                  fine(i,jj+2,k,n), fine(i,jj+3,k,n),
                                  fine(i,jj+4,k,n), fine(i,jj+5,k,n),
                                  1, true, true, crse(i,j,k,n));
        } else if (j == big_end) {
            WENO_Restrict_Unified(fine(i,jj+3,k,n),   fine(i,jj+2,k,n),
                                  fine(i,jj+1,k,n), fine(i,jj,k,n),
                                  fine(i,jj-1,k,n), fine(i,jj-2,k,n),
                                  1, true, true, crse(i,j,k,n));
        }
    } else {
        amrex::Abort("weno_restrict_y: how did this happen?");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_z (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine, IntVect const& ratio) noexcept
{
    int small_end = crse.begin.z, big_end = crse.end.z-1;
    int kk = k*ratio[2];
    if (ratio[2] == 2) {
        if (k > small_end && k < big_end) {
            WENO_Restrict_Unified(fine(i,j,kk-2,n), fine(i,j,kk-1,n),
                                  fine(i,j,kk,n),   fine(i,j,kk+1,n),
                                  fine(i,j,kk+2,n), fine(i,j,kk+3,n),
                                  0, false, false, crse(i,j,k,n));
        } else if (k == small_end) {
            WENO_Restrict_Unified(fine(i,j,kk,n),   fine(i,j,kk+1,n),
                                  fine(i,j,kk+2,n), fine(i,j,kk+3,n),
                                  fine(i,j,kk+4,n), fine(i,j,kk+5,n),
                                  0, true, true, crse(i,j,k,n));
        } else if (k == big_end) {
            WENO_Restrict_Unified(fine(i,j,kk+1,n), fine(i,j,kk,n),
                                  fine(i,j,kk-1,n), fine(i,j,kk-2,n),
                                  fine(i,j,kk-3,n), fine(i,j,kk-4,n),
                                  0, true, true, crse(i,j,k,n));
        }
    } else if (ratio[2] == 4) {
        if (k > small_end && k < big_end) {
            WENO_Restrict_Unified(fine(i,j,kk-1,n), fine(i,j,kk,n),
                                  fine(i,j,kk+1,n), fine(i,j,kk+2,n),
                                  fine(i,j,kk+3,n), fine(i,j,kk+4,n),
                                  1, false, false, crse(i,j,k,n));
        } else if (k == small_end) {
            WENO_Restrict_Unified(fine(i,j,kk,n),   fine(i,j,kk+1,n),
                                  fine(i,j,kk+2,n), fine(i,j,kk+3,n),
                                  fine(i,j,kk+4,n), fine(i,j,kk+5,n),
                                  1, true, true, crse(i,j,k,n));
        } else if (k == big_end) {
            WENO_Restrict_Unified(fine(i,j,kk+3,n),   fine(i,j,kk+2,n),
                                  fine(i,j,kk+1,n), fine(i,j,kk,n),
                                  fine(i,j,kk-1,n), fine(i,j,kk-2,n),
                                  1, true, true, crse(i,j,k,n));
        }
    } else {
        amrex::Abort("weno_restrict_z: how did this happen?");
    }
}

}
#endif
