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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void WENO_Interp_JS6(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
                       const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
                       amrex::Real& v_edge)
{
    using namespace amrex;
   
    // WENO-JS linear weights
    constexpr Real d0 = 3.0/256.0;
    constexpr Real d1 = 99.0/256.0;
    constexpr Real d2 = 693.0/1280.0;
    constexpr Real d3 = 77.0/1280.0;

    // Compute smoothness indicators
    Real beta0 = 1.0/3.0*(
                    29.0*v0*v2 - 49.0*v0*v1 - 73.0*v1*v2 + 
                    10.0*std::pow(v0,2) + 61.0*std::pow(v1,2) + 22.0*std::pow(v2,2)
    );
    Real beta1 = 1.0/3.0*(
                    11.0*v1*v3 - 19.0*v1*v2 - 31.0*v2*v3 + 
                    4.0*std::pow(v1,2) + 25.0*std::pow(v2,2) + 10.0*std::pow(v3,2)
    );
    Real beta2 = 1.0/3.0*(
                    5.0*v2*v4 - 13.0*v2*v3 - 13.0*v3*v4 + 
                    4.0*std::pow(v2,2) + 13.0*std::pow(v3,2) + 4.0*std::pow(v4,2)
    );
    Real beta3 = 1.0/3.0*(
                    11.0*v3*v5 - 31.0*v3*v4 - 19.0*v4*v5 + 
                    10.0*std::pow(v3,2) + 25.0*std::pow(v4,2) + 4.0*std::pow(v5,2)
    );

    // small value for computing non-linear weights
    constexpr Real eps = 1.0e-6;

    // alpha values WENO-JS
    Real alpha0 = d0/std::pow(eps+beta0, 2);
    Real alpha1 = d1/std::pow(eps+beta1, 2);
    Real alpha2 = d2/std::pow(eps+beta2, 2);
    Real alpha3 = d3/std::pow(eps+beta3, 2);

    // Normalize weights
    Real sum_alpha = alpha0 + alpha1 + alpha2 + alpha3;
    Real omega0 = alpha0/sum_alpha;
    Real omega1 = alpha1/sum_alpha;
    Real omega2 = alpha2/sum_alpha;
    Real omega3 = alpha3/sum_alpha;

    // Candidate stencils
    Real q0 = (1.0/32.0)*(21.0*v0 - 66.0*v1 + 77.0*v2);
    Real q1 = (1.0/32.0)*(-3.0*v1 + 14.0*v2 + 21.0*v3);
    Real q2 = (1.0/32.0)*( 5.0*v2 + 30.0*v3  - 3.0*v4);
    Real q3 = (1.0/32.0)*(45.0*v3 - 18.0*v4  + 5.0*v5);

    // WENO-JS reconstruction
    v_edge = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_x (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    // index 
    int ii = amrex::coarsen(i,2);
    int s = 2 * (i - ii*2) - 1; // i=ii*2, s = -1; i=ii*2+1, s = 1

    WENO_Interp_JS6(crse(ii+s*3,j,k,n), crse(ii+s*2,j,k,n), 
                    crse(ii+s*1,j,k,n), crse(ii,j,k,n),
                    crse(ii-s*1,j,k,n), crse(ii-s*2,j,k,n),
                    fine(i,j,k,n));

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_y (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    // index 
    int jj = amrex::coarsen(j,2);
    int s = 2 * (j - jj*2) - 1; // j=jj*2, s = -1; j=jj*2+1, s = 1

    WENO_Interp_JS6(crse(i,jj+s*3,k,n), crse(i,jj+s*2,k,n), 
                    crse(i,jj+s*1,k,n), crse(i,jj,k,n),
                    crse(i,jj-s*1,k,n), crse(i,jj-s*2,k,n),
                    fine(i,j,k,n));

}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void central_weno_interp_z (int i, int j, int k, int n, 
                            Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    // index 
    int kk = amrex::coarsen(i,2);
    int s = 2 * (k - kk*2) - 1; // k=kk*2, s = -1; k=kk*2+1, s = 1

    WENO_Interp_JS6(crse(i,j,kk+s*3,n), crse(i,j,kk+s*2,n), 
                    crse(i,j,kk+s*1,n), crse(i,j,kk,n),
                    crse(i,j,kk-s*1,n), crse(i,j,kk-s*2,n),
                    fine(i,j,k,n));

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void WENO_Restrict_JS6_c(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
                       const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
                       amrex::Real& v_edge)
{
    using namespace amrex;

    // WENO-JS linear weights
    constexpr Real d0 =  1.0/32.0;
    constexpr Real d1 = 15.0/32.0;
    constexpr Real d2 = 15.0/32.0;
    constexpr Real d3 =  1.0/32.0;

    // Compute smoothness indicators
    Real beta0 = 1.0/12.0*(
                    74.0*v0*v2 - 124.0*v0*v1 - 196.0*v1*v2 + 
                    25.0*std::pow(v0,2) + 160.0*std::pow(v1,2) + 61.0*std::pow(v2,2)
    );
    Real beta1 = 1.0/12.0*(
                    26.0*v1*v3 - 52.0*v1*v2 - 76.0*v2*v3 + 
                    13.0*std::pow(v1,2) + 64.0*std::pow(v2,2) + 25.0*std::pow(v3,2)
    );
    Real beta2 = 1.0/12.0*(
                    26.0*v2*v4 - 76.0*v2*v3 - 52.0*v3*v4 +
                    25.0*std::pow(v2,2) + 64.0*std::pow(v3,2) + 13.0*std::pow(v4,2)
    );
    Real beta3 = 1.0/12.0*(
                    74.0*v3*v5 - 196.0*v3*v4 - 124.0*v4*v5 + 
                    61.0*std::pow(v3,2) + 160.0*std::pow(v4,2) + 25.0*std::pow(v5,2)
    );

    // small value for computing non-linear weights
    constexpr Real eps = 1.0e-6;

    // alpha values WENO-JS
    Real alpha0 = d0/std::pow(eps+beta0, 2);
    Real alpha1 = d1/std::pow(eps+beta1, 2);
    Real alpha2 = d2/std::pow(eps+beta2, 2);
    Real alpha3 = d3/std::pow(eps+beta2, 2);

    // Normalize weights
    Real sum_alpha = alpha0 + alpha1 + alpha2 + alpha3;
    Real omega0 = alpha0/sum_alpha;
    Real omega1 = alpha1/sum_alpha;
    Real omega2 = alpha2/sum_alpha;
    Real omega3 = alpha3/sum_alpha;

    // Candidate stencils
    Real q0 = (1.0/8.0)*( 3.0*v0 - 10.0*v1 + 15.0*v2);
    Real q1 = (1.0/8.0)*(-1.0*v1 +  6.0*v2 + 3.0*v3);
    Real q2 = (1.0/8.0)*( 3.0*v2 +  6.0*v3 - 1.0*v4);
    Real q3 = (1.0/8.0)*(15.0*v3 - 10.0*v4 + 3.0*v5);

    // WENO-JS reconstruction
    v_edge = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3;

}

// 6-point upwind WENO is very unstable
// AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
// void WENO_Restrict_JS6_b(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
//                        const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
//                        amrex::Real& v_edge)
// {
//     using namespace amrex;

//     // WENO-JS linear weights
//     constexpr Real d0 = 21.0/32.0;
//     constexpr Real d1 = 63.0/160.0;
//     constexpr Real d2 = -9.0/160.0;
//     constexpr Real d3 =  1.0/160.0;

//     // Compute smoothness indicators
//     Real beta0 = 1.0/12.0*(
//                     26.0*v0*v2 - 76.0*v0*v1 - 52.0*v1*v2 + 
//                     25.0*std::pow(v0,2) + 64.0*std::pow(v1,2) + 13.0*std::pow(v2,2)
//     );
//     Real beta1 = 1.0/12.0*(
//                     74.0*v1*v3 - 196.0*v1*v2 - 124.0*v2*v3 + 
//                     61.0*std::pow(v1,2) + 160.0*std::pow(v2,2) + 25.0*std::pow(v3,2)
//     );
//     Real beta2 = 1.0/12.0*(
//                     170.0*v2*v4 - 412.0*v2*v3 - 292.0*v3*v4 +
//                     121.0*std::pow(v2,2) + 352.0*std::pow(v3,2) + 61.0*std::pow(v4,2)
//     );
//     Real beta3 = 1.0/12.0*(
//                     314.0*v3*v5 - 724.0*v3*v4 - 556.0*v4*v5 + 
//                     205.0*std::pow(v3,2) + 640.0*std::pow(v4,2) + 121.0*std::pow(v5,2)
//     );

//     // small value for computing non-linear weights
//     constexpr Real eps = 1.0e-6;

//     // alpha values WENO-JS
//     Real alpha0 = d0/std::pow(eps+beta0, 2);
//     Real alpha1 = d1/std::pow(eps+beta1, 2);
//     Real alpha2 = d2/std::pow(eps+beta2, 2);
//     Real alpha3 = d3/std::pow(eps+beta2, 2);

//     // Normalize weights
//     Real sum_alpha = alpha0 + alpha1 + alpha2 + alpha3;
//     Real omega0 = alpha0/sum_alpha;
//     Real omega1 = alpha1/sum_alpha;
//     Real omega2 = alpha2/sum_alpha;
//     Real omega3 = alpha3/sum_alpha;

//     // Candidate stencils
//     Real q0 = (1.0/8.0)*( 3.0*v0 + 6.0*v1 - 1.0*v2);
//     Real q1 = (1.0/8.0)*(15.0*v1 -10.0*v2 + 3.0*v3);
//     Real q2 = (1.0/8.0)*(35.0*v2 -42.0*v3 + 15.0*v4);
//     Real q3 = (1.0/8.0)*(63.0*v3 -90.0*v4 + 35.0*v5);

//     // WENO-JS reconstruction
//     v_edge = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3;
// }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void WENO_Restrict_JS6_b(const amrex::Real& v0, const amrex::Real& v1, const amrex::Real& v2, 
                       const amrex::Real& v3, const amrex::Real& v4, const amrex::Real& v5, 
                       amrex::Real& v_edge)
{
    using namespace amrex;

    // 2nd-order Lagragnge interpolation
    v_edge = (1.0/8.0)*( 3.0*v0 + 6.0*v1 - 1.0*v2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_cons_restrict(const amrex::Real& v0, const amrex::Real& v1, amrex::Real& v_edge)
{
    v_edge = (v0+v1)/2.0;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_x (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine) noexcept
{
    int small_end = crse.begin.x, big_end = crse.end.x-1;
    int ii = i*2;
    if (i > small_end && i < big_end) {
        WENO_Restrict_JS6_c(fine(ii-2,j,k,n), fine(ii-1,j,k,n),
                            fine(ii,j,k,n),   fine(ii+1,j,k,n),
                            fine(ii+2,j,k,n), fine(ii+3,j,k,n),
                            crse(i,j,k,n));
    } else if (i == small_end) {
        WENO_Restrict_JS6_b(fine(ii,j,k,n),   fine(ii+1,j,k,n),
                            fine(ii+2,j,k,n), fine(ii+3,j,k,n),
                            fine(ii+4,j,k,n), fine(ii+5,j,k,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(ii,j,k,n), fine(ii+1,j,k,n), crse(i,j,k,n));
    } else if (i == big_end) {
        WENO_Restrict_JS6_b(fine(ii+1,j,k,n), fine(ii,j,k,n),
                            fine(ii-1,j,k,n), fine(ii-2,j,k,n),
                            fine(ii-3,j,k,n), fine(ii-4,j,k,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(ii,j,k,n), fine(ii+1,j,k,n), crse(i,j,k,n));
    } else {
        amrex::Abort("weno_restrict_x: how did this happen?");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_y (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine) noexcept
{
    int small_end = crse.begin.y, big_end = crse.end.y-1;
    int jj = j*2;
    if (j > small_end && j < big_end) {
        WENO_Restrict_JS6_c(fine(i,jj-2,k,n), fine(i,jj-1,k,n),
                            fine(i,jj,k,n),   fine(i,jj+1,k,n),
                            fine(i,jj+2,k,n), fine(i,jj+3,k,n),
                            crse(i,j,k,n));
    } else if (j == small_end) {
        WENO_Restrict_JS6_b(fine(i,jj,k,n),   fine(i,jj+1,k,n),
                            fine(i,jj+2,k,n), fine(i,jj+3,k,n),
                            fine(i,jj+4,k,n), fine(i,jj+5,k,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(i,jj,k,n), fine(i,jj+1,k,n), crse(i,j,k,n));
    } else if (j == big_end) {
        WENO_Restrict_JS6_b(fine(i,jj+1,k,n), fine(i,jj,k,n),
                            fine(i,jj-1,k,n), fine(i,jj-2,k,n),
                            fine(i,jj-3,k,n), fine(i,jj-4,k,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(i,jj,k,n), fine(i,jj+1,k,n), crse(i,j,k,n));
    } else {
        amrex::Abort("weno_restrict_y:  how did this happen?");
    }
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void weno_restrict_z (int i, int j, int k, int n, Array4<Real> const& crse,
                      Array4<Real const> const& fine) noexcept
{
    int small_end = crse.begin.z, big_end = crse.end.z-1;
    int kk = k*2;
    if (k > small_end && k < big_end) {
        WENO_Restrict_JS6_c(fine(i,j,kk-2,n), fine(i,j,kk-1,n),
                            fine(i,j,kk,n),   fine(i,j,kk+1,n),
                            fine(i,j,kk+2,n), fine(i,j,kk+3,n),
                            crse(i,j,k,n));
    } else if (k == small_end) {
        WENO_Restrict_JS6_b(fine(i,j,kk,n),   fine(i,j,kk+1,n),
                            fine(i,j,kk+2,n), fine(i,j,kk+3,n),
                            fine(i,j,kk+4,n), fine(i,j,kk+5,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(i,j,kk,n), fine(i,j,kk+1,n), crse(i,j,k,n));
    } else if (k == big_end) {
        WENO_Restrict_JS6_b(fine(i,j,kk+1,n), fine(i,j,kk,n),
                            fine(i,j,kk-1,n), fine(i,j,kk-2,n),
                            fine(i,j,kk-3,n), fine(i,j,kk-4,n),
                            crse(i,j,k,n));
        // cell_cons_restrict(fine(i,j,kk,n), fine(i,j,kk+1,n), crse(i,j,k,n));
    } else {
        amrex::Abort("weno_restrict_z:  how did this happen?");
    }
}

}
#endif
