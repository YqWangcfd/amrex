
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_ParmParse.H>

#include <climits>
#include <cmath>
#include <IndexMacro.H>

namespace amrex {

/*
 * PCInterp, NodeBilinear, FaceLinear, CellConservativeLinear and
 * CellBilinear are supported for all dimensions on cpu and gpu.
 *
 * CellConservativeProtected only works in 2D and 3D on cpu and gpu
 * and assumes that ratio > 1 in all directions
 *
 * CellQuadratic only works in 2D and 3D on cpu and gpu.
 *
 * CellQuartic works in 1D, 2D and 3D on cpu and gpu with ref ratio of 2
 *
 * CellConservativeQuartic only works with ref ratio of 2 on cpu and gpu.
 *
 * FaceConservativeLinear works in 2D and 3D on cpu and gpu.
 *
 * FaceDivFree works in 2D and 3D on cpu and gpu.
 * The algorithm is restricted to ref ratio of 2.
 */

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
PCInterp                  pc_interp;
NodeBilinear              node_bilinear_interp;
FaceLinear                face_linear_interp;
FaceConservativeLinear    face_cons_linear_interp;
FaceDivFree               face_divfree_interp;
CellConservativeLinear    lincc_interp;
CellConservativeLinear    cell_cons_interp(false);
CellConservativeProtected protected_interp;
CellConservativeQuartic   quartic_interp;
CellBilinear              cell_bilinear_interp;
CellQuadratic             quadratic_interp;
CellQuartic               cell_quartic_interp;
CellWENO                  cell_weno_interp;
Mortar2D                  mortar_interp_orderRef;
Mortar2D                  mortar_interp_scaleRef(Mortar2D::Type::ScaleRef);
HermiteWENO2D             hermite_weno_interp;

namespace{

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real hweno_pressure_ideal (Real rho, Real mx, Real my, Real mz, Real e, Real gamma) noexcept
    {
        const Real rho_safe = amrex::max(rho, Real(1.0e-10));
        const Real kinetic = Real(0.5) * (mx*mx + my*my + mz*mz) / rho_safe;
        return (gamma - Real(1.0)) * (e - kinetic);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void hweno_apply_positivity_cell (Array4<Real> const& arr, int i, int j, int k, int nvar,
                                      Real eos_gamma, Real eps_rho, Real eps_p) noexcept
    {
        if (nvar <= ME) { return; }

        auto sd_weight = [] AMREX_GPU_HOST_DEVICE (int off) noexcept -> Real
        {
            const int ix = off % sd_ORDER;
            const int iy = (off / sd_ORDER) % sd_ORDER;
#if (AMREX_SPACEDIM == 3)
            const int iz = off / (sd_ORDER * sd_ORDER);
            return SDwgh1D[ix] * SDwgh1D[iy] * SDwgh1D[iz];
#else
            return SDwgh1D[ix] * SDwgh1D[iy];
#endif
        };

        Real rho_bar = Real(0.0);
        Real mx_bar  = Real(0.0);
        Real my_bar  = Real(0.0);
        Real mz_bar  = Real(0.0);
        Real e_bar   = Real(0.0);
        for (int off = 0; off < HermiteWENO2D::sd_space_hweno; ++off) {
            const Real w = sd_weight(off);
            rho_bar += w * arr(i,j,k,MRHO*HermiteWENO2D::sd_space_hweno + off);
            mx_bar  += w * arr(i,j,k,MU  *HermiteWENO2D::sd_space_hweno + off);
            my_bar  += w * arr(i,j,k,MV  *HermiteWENO2D::sd_space_hweno + off);
            mz_bar  += w * arr(i,j,k,MW  *HermiteWENO2D::sd_space_hweno + off);
            e_bar   += w * arr(i,j,k,ME  *HermiteWENO2D::sd_space_hweno + off);
        }

        Real rho_min = AMREX_REAL_MAX;
        for (int off = 0; off < HermiteWENO2D::sd_space_hweno; ++off) {
            rho_min = amrex::min(rho_min, arr(i,j,k,MRHO*HermiteWENO2D::sd_space_hweno + off));
        }

        Real theta1 = Real(1.0);
        if (rho_min < eps_rho) {
            const Real denom = rho_bar - rho_min;
            if (denom > Real(1.0e-30)) {
                theta1 = amrex::min(Real(1.0), (rho_bar - eps_rho) / denom);
            } else {
                theta1 = Real(0.0);
            }
        }

        Real theta2 = Real(1.0);
        for (int off = 0; off < HermiteWENO2D::sd_space_hweno; ++off) {
            const Real rho_raw = arr(i,j,k,MRHO*HermiteWENO2D::sd_space_hweno + off);
            const Real mx_raw  = arr(i,j,k,MU  *HermiteWENO2D::sd_space_hweno + off);
            const Real my_raw  = arr(i,j,k,MV  *HermiteWENO2D::sd_space_hweno + off);
            const Real mz_raw  = arr(i,j,k,MW  *HermiteWENO2D::sd_space_hweno + off);
            const Real e_raw   = arr(i,j,k,ME  *HermiteWENO2D::sd_space_hweno + off);

            const Real rho_hat = rho_bar + theta1 * (rho_raw - rho_bar);
            const Real p_hat = hweno_pressure_ideal(rho_hat, mx_raw, my_raw, mz_raw, e_raw, eos_gamma);

            if (!std::isfinite(p_hat) || rho_hat < eps_rho || p_hat < eps_p) {
                Real lo = Real(0.0);
                Real hi = Real(1.0);
                for (int iter = 0; iter < 30; ++iter) {
                    const Real mid = Real(0.5) * (lo + hi);
                    const Real rho_mid = rho_bar + mid * (rho_hat - rho_bar);
                    const Real mx_mid  = mx_bar  + mid * (mx_raw  - mx_bar);
                    const Real my_mid  = my_bar  + mid * (my_raw  - my_bar);
                    const Real mz_mid  = mz_bar  + mid * (mz_raw  - mz_bar);
                    const Real e_mid   = e_bar   + mid * (e_raw   - e_bar);
                    const Real p_mid = hweno_pressure_ideal(rho_mid, mx_mid, my_mid, mz_mid, e_mid, eos_gamma);
                    if (rho_mid >= eps_rho && std::isfinite(p_mid) && p_mid >= eps_p) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                theta2 = amrex::min(theta2, lo);
            }
        }

        if (theta1 < Real(1.0) || theta2 < Real(1.0)) {
            for (int n = 0; n < nvar; ++n) {
                Real ubar = Real(0.0);
                for (int off = 0; off < HermiteWENO2D::sd_space_hweno; ++off) {
                    ubar += sd_weight(off) * arr(i,j,k,n*HermiteWENO2D::sd_space_hweno + off);
                }
                for (int off = 0; off < HermiteWENO2D::sd_space_hweno; ++off) {
                    const Real raw = arr(i,j,k,n*HermiteWENO2D::sd_space_hweno + off);
                    const Real uhat = (n == MRHO) ? (ubar + theta1 * (raw - ubar)) : raw;
                    arr(i,j,k,n*HermiteWENO2D::sd_space_hweno + off) = ubar + theta2 * (uhat - ubar);
                }
            }
        }
    }
}

Box
NodeBilinear::CoarseBox (const Box& fine,
                         int        ratio)
{
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

Box
NodeBilinear::CoarseBox (const Box&     fine,
                         const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

void
NodeBilinear::interp (const FArrayBox&  crse,
                      int               crse_comp,
                      FArrayBox&        fine,
                      int               fine_comp,
                      int               ncomp,
                      const Box&        fine_region,
                      const IntVect&    ratio,
                      const Geometry& /*crse_geom */,
                      const Geometry& /*fine_geom */,
                      Vector<BCRec> const& /*bcr*/,
                      int               /*actual_comp*/,
                      int               /*actual_state*/,
                      RunOn             runon)
{
    BL_PROFILE("NodeBilinear::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
    {
        mf_nodebilin_interp(i,j,k,n, finearr, fine_comp, crsearr, crse_comp, ratio);
    });
}

Box
FaceLinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

Box
FaceLinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (b.type(i) == IndexType::NODE && b.length(i) < 2) {
            // Don't want degenerate boxes in nodal direction.
            b.growHi(i,1);
        }
    }
    return b;
}

void
FaceLinear::interp (const FArrayBox&  crse,
                    int               crse_comp,
                    FArrayBox&        fine,
                    int               fine_comp,
                    int               ncomp,
                    const Box&        fine_region,
                    const IntVect&    ratio,
                    const Geometry& /* crse_geom */,
                    const Geometry& /* fine_geom */,
                    Vector<BCRec> const& /*bcr*/,
                    int              /* actual_comp*/,
                    int               /*actual_state*/,
                    RunOn             runon)
{
    //
    // This version is called from FillPatchInterp which is called by
    //      InterpFromCoarseLevel in AMReX_FillPatchUtil_I.H
    //
    // It assumes no existing fine values that need to be preserved (unlike interp_face below)
    //
    // Inside each call to face_linear_interp_* (in AMRex_Interp_*D_C.H), we do:
    //  * on fine faces which overlie crse faces, the fine value is set to the crse value (piecewise constant)
    //  * on fine faces which are between two crse faces, the fine value is set to the average of the crse values (linear)
    //
    BL_PROFILE("FaceLinear::interp()");

    AMREX_ASSERT(AMREX_D_TERM(fine_region.type(0),+fine_region.type(1),+fine_region.type(2)) == 1);

    Array4<Real> const& fine_arr = fine.array(fine_comp);
    Array4<Real const> const& crse_arr = crse.const_array(crse_comp);

    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_x(i,j,k,n,fine_arr,crse_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_y(i,j,k,n,fine_arr,crse_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_z(i,j,k,n,fine_arr,crse_arr,ratio);
        });
    }
#endif
#endif
}

void
FaceLinear::interp_face (const FArrayBox&  crse,
                         const int         crse_comp,
                         FArrayBox&        fine,
                         const int         fine_comp,
                         const int         ncomp,
                         const Box&        fine_region,
                         const IntVect&    ratio,
                         const IArrayBox&  solve_mask,
                         const Geometry& /*crse_geom */,
                         const Geometry& /*fine_geom */,
                         Vector<BCRec> const& /*bcr*/,
                         const int         /*bccomp*/,
                         RunOn             runon)
{
    //
    // This version is called from InterpFace which is called from the version FillPatchTwoLevels_doit
    //      that takes a single MF (in AMReX_FillPatchUtil_I.H)
    //
    // It assumes there are existing fine values which we want to preserve (unlike interp above)
    //
    // We do the interpolation in two steps:
    //   1) face_linear_face_interp_*: on fine faces which overlie crse faces, the fine value is set to the crse value (piecewise constant) ONLY IF
    //      there is not already fine data there
    //   2) face_linear_interp_*: on fine faces which are between two crse faces, the fine value is set to the average of the values
    //      on the faces overlying -- this uses only the results of step 1, it does not take the crse values
    //
    BL_PROFILE("FaceLinear::interp_face()");

    AMREX_ASSERT(AMREX_D_TERM(fine_region.type(0),+fine_region.type(1),+fine_region.type(2)) == 1);

    Array4<Real> const& fine_arr = fine.array(fine_comp);
    Array4<Real const> const& crse_arr = crse.const_array(crse_comp);
    Array4<const int> mask_arr;
    if (solve_mask.isAllocated()) {
        mask_arr = solve_mask.const_array();
    }

    //
    // Fill fine ghost faces with piecewise-constant interpolation of coarse data.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_x(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_y(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_z(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#endif
#endif

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_x(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_y(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_z(i,j,k,n,fine_arr,ratio);
        });
    }
#endif
#endif
}

void FaceLinear::interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                             const int         crse_comp,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                             const int         fine_comp,
                             const int         ncomp,
                             const Box&        fine_region,
                             const IntVect&    ratio,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                             const Geometry&   /*crse_geom*/,
                             const Geometry&   /*fine_geom*/,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                             const int         /*actual_comp*/,
                             const int         /*actual_state*/,
                             const RunOn       runon)
{
    //
    // This version is called from FillPatchTwoLevels_doit (that takes an Array of MF*) in AMReX_FillPatchUtil_I.H
    //
    // It assumes there are existing fine values which we want to preserve (like face_interp, unlike interp above)
    //
    // We do the interpolation in two steps:
    //   1) face_linear_face_interp_*: on fine faces which overlie crse faces, the fine value is set to the crse value (piecewise constant) ONLY IF
    //      there is not already fine data there
    //   2) face_linear_interp_*: on fine faces which are between two crse faces, the fine value is set to the average of the values
    //      on the faces overlying -- this uses only the results of step 1, it does not take the crse values
    //
    BL_PROFILE("FaceLinear::interp_arr()");

    Array<IndexType, AMREX_SPACEDIM> types;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        { types[d].set(d); }

    GpuArray<Array4<const Real>, AMREX_SPACEDIM> crse_arr;
    GpuArray<Array4<Real>, AMREX_SPACEDIM> fine_arr;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> mask_arr;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        crse_arr[d] = crse[d]->const_array(crse_comp);
        fine_arr[d] = fine[d]->array(fine_comp);
        if (solve_mask[d] != nullptr)
            { mask_arr[d] = solve_mask[d]->const_array(0); }
    }

    //
    // Fill fine ghost faces with piecewise-constant interpolation of coarse data.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    // Fuse the launches, 1 for each dimension, into a single launch.
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_x(i,j,k,n,fine_arr[0],crse_arr[0],mask_arr[0],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_y(i,j,k,n,fine_arr[1],crse_arr[1],mask_arr[1],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_z(i,j,k,n,fine_arr[2],crse_arr[2],mask_arr[2],ratio);
                      }
                  });
              });

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_x(i,j,k,n,fine_arr[0],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_y(i,j,k,n,fine_arr[1],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_z(i,j,k,n,fine_arr[2],ratio);
                      }
                  });
              });
}

Box
FaceConservativeLinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

Box
FaceConservativeLinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    IntVect ng(1);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if ( (fine.type(i) == IndexType::NODE) || (ratio[i] == 1) ) {
            ng[i] = 0;
        }
    }
    Box b = amrex::coarsen(fine,ratio); b.grow(ng);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (b.type(i) == IndexType::NODE) {
            if (b.type(i) == IndexType::NODE && b.length(i) < 2) {
                // Don't want degenerate boxes in nodal direction.
                b.growHi(i,1);
            }
        }
    }
    return b;
}

void
FaceConservativeLinear::interp (const FArrayBox&     crse,
                                int                  crse_comp,
                                FArrayBox&           fine,
                                int                  fine_comp,
                                int                  ncomp,
                                const Box&           fine_region,
                                const IntVect&       ratio,
                                const Geometry&      crse_geom,
                                const Geometry&      fine_geom,
                                Vector<BCRec> const& bcr,
                                int                  /*actual_comp*/,
                                int                  /*actual_state*/,
                                RunOn                runon)
{
    //
    // This version is called from FillPatchInterp which is called by
    //      InterpFromCoarseLevel in AMReX_FillPatchUtil_I.H
    //
    // It assumes no existing fine values that need to be preserved thus does not send a mask to interp_face
    //
    BL_PROFILE("FaceConservativeLinear::interp()");

    AMREX_ASSERT(AMREX_D_TERM(fine_region.type(0),+fine_region.type(1),+fine_region.type(2)) == 1);

    // We intentionally do not allocate the mask so that all faces are filled from coarse values
    IArrayBox dummy_mask;
    int bccomp = 0; // This is also a dummy -- it's not used
    interp_face(crse,crse_comp,fine,fine_comp,ncomp,fine_region,ratio,dummy_mask,
                crse_geom,fine_geom,bcr,bccomp,runon);
}

void
FaceConservativeLinear::interp_face (const FArrayBox&       crse,
                                     const int              crse_comp,
                                     FArrayBox&             fine,
                                     const int              fine_comp,
                                     const int              ncomp,
                                     const Box&             fine_region,
                                     const IntVect&         ratio,
                                     const IArrayBox&       solve_mask,
                                     const Geometry&        crse_geom,
                                     const Geometry&      /*fine_geom */,
                                     Vector<BCRec> const& /*bcr*/,
                                     const int            /*bccomp*/,
                                     RunOn                  runon)
{
    //
    // This version is called from InterpFace which is called from the version FillPatchTwoLevels_doit
    //      that takes a single MF (in AMReX_FillPatchUtil_I.H)
    //
    // It assumes there are existing fine values which we want to preserve (unlike interp above)
    //
    // We do the interpolation in two steps:
    //   1) face_cons_linear_face_interp: on fine faces which overlie crse faces, slopes are computed (linear in 2d, bilinear in 3d)
    //      and the fine value is over-written ONLY IF there is not already fine data there (assuming the mask is used)
    //   2) face_linear_interp_*: on fine faces which are between two crse faces, the fine value is set to the average of the values
    //      on the faces overlying -- this uses only the results of step 1
    //      NOTE: we use the same routines as used by FaceLinear since this interpolation is only in the normal direction
    //
    BL_PROFILE("FaceConservativeLinear::interp_face()");

    AMREX_ASSERT(AMREX_D_TERM(fine_region.type(0),+fine_region.type(1),+fine_region.type(2)) == 1);
    Array4<Real> const& fine_arr = fine.array(fine_comp);
    Array4<Real const> const& crse_arr = crse.const_array(crse_comp);
    Array4<const int> mask_arr;
    if (solve_mask.isAllocated()) {
        mask_arr = solve_mask.const_array();
    }

    // We don't need to worry about face-based domain because this is only used in the tangential interpolation
    Box per_grown_domain = crse_geom.Domain();
    for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
        if (crse_geom.isPeriodic(dim)) {
            per_grown_domain.grow(dim,1);
        }
    }

    //
    // Fill fine ghost faces with interpolation of coarse data that is conservative linear
    //      in the tangential direction.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_cons_linear_face_interp(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio,per_grown_domain,0);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_cons_linear_face_interp(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio,per_grown_domain,1);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_cons_linear_face_interp(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio,per_grown_domain,2);
        });
    }
#endif
#endif

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_x(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_y(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_z(i,j,k,n,fine_arr,ratio);
        });
    }
#endif
#endif
}

void FaceConservativeLinear::interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                                         const int         crse_comp,
                                         Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                                         const int         fine_comp,
                                         const int         ncomp,
                                         const Box&        fine_region,
                                         const IntVect&    ratio,
                                         Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                                         const Geometry&   crse_geom,
                                         const Geometry&   /*fine_geom*/,
                                         Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                                         const int         /*actual_comp*/,
                                         const int         /*actual_state*/,
                                         const RunOn       runon)
{
    //
    // This version is called from FillPatchTwoLevels_doit (that takes an Array of MF*) in AMReX_FillPatchUtil_I.H
    //
    // It assumes there are existing fine values which we want to preserve (like face_interp, unlike interp above)
    //
    // We do the interpolation in two steps:
    //   1) face_cons_linear_face_interp_*: on fine faces which overlie crse faces, we compute tangential slopes
    //      to compute the fine values (linear in 2d, bilinear in 3d) ONLY IF there is not already fine data there
    //   2) face_cons_linear_interp_*: on fine faces which are between two crse faces, the fine value is set to the average of the values
    //      on the faces overlying -- this uses only the results of step 1, it does not take the crse values
    //      NOTE: here we use the same routines as used by FaceLinear since this interpolation is only in the normal direction
    //
    BL_PROFILE("FaceConservativeLinear::interp_arr()");

    Array<IndexType, AMREX_SPACEDIM> types;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        { types[d].set(d); }

    GpuArray<Array4<const Real>, AMREX_SPACEDIM> crse_arr;
    GpuArray<Array4<Real>, AMREX_SPACEDIM> fine_arr;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> mask_arr;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        crse_arr[d] = crse[d]->const_array(crse_comp);
        fine_arr[d] = fine[d]->array(fine_comp);
        if (solve_mask[d] != nullptr)
            { mask_arr[d] = solve_mask[d]->const_array(0); }
    }

    // We don't need to worry about face-based domain because this is only used in the tangential interpolation
    Box per_grown_domain = crse_geom.Domain();
    for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
        if (crse_geom.isPeriodic(dim)) {
            per_grown_domain.grow(dim,1);
        }
    }

    //
    // Fill fine ghost faces with interpolation of coarse data that is conservative linear
    //      in the tangential direction.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    // Fuse the launches, 1 for each dimension, into a single launch.
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_cons_linear_face_interp(i,j,k,n,fine_arr[0],crse_arr[0],mask_arr[0],ratio,per_grown_domain,0);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_cons_linear_face_interp(i,j,k,n,fine_arr[1],crse_arr[1],mask_arr[1],ratio,per_grown_domain,1);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_cons_linear_face_interp(i,j,k,n,fine_arr[2],crse_arr[2],mask_arr[2],ratio,per_grown_domain,2);
                      }
                  });
              });

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_x(i,j,k,n,fine_arr[0],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_y(i,j,k,n,fine_arr[1],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_z(i,j,k,n,fine_arr[2],ratio);
                      }
                  });
              });
}

Box
CellBilinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

Box
CellBilinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    const int* lo = fine.loVect();
    const int* hi = fine.hiVect();

    Box crse(amrex::coarsen(fine,ratio));
    const int* clo = crse.loVect();
    const int* chi = crse.hiVect();

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if ((lo[i]-clo[i]*ratio[i])*2 < ratio[i]) {
            crse.growLo(i,1);
        }
        if ((hi[i]-chi[i]*ratio[i])*2 >= ratio[i]) {
            crse.growHi(i,1);
        }
    }
    return crse;
}

void
CellBilinear::interp (const FArrayBox&  crsefab,
                      int               crse_comp,
                      FArrayBox&        finefab,
                      int               fine_comp,
                      int               ncomp,
                      const Box&        fine_region,
                      const IntVect &   ratio,
                      const Geometry& /*crse_geom*/,
                      const Geometry& /*fine_geom*/,
                      Vector<BCRec> const& /*bcr*/,
                      int               /*actual_comp*/,
                      int               /*actual_state*/,
                      RunOn             runon)
{
    BL_PROFILE("CellBilinear::interp()");

    auto const& crse = crsefab.const_array();
    auto const& fine = finefab.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
    {
        mf_cell_bilin_interp(i,j,k,n, fine, fine_comp, crse, crse_comp, ratio);
    });
}


CellConservativeLinear::CellConservativeLinear (bool do_linear_limiting_)
    : do_linear_limiting(do_linear_limiting_)
{}

Box
CellConservativeLinear::CoarseBox (const Box&     fine,
                                   const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
        if (ratio[dim] > 1) {
            crse.grow(dim,1);
        }
    }
    return crse;
}

Box
CellConservativeLinear::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
CellConservativeLinear::interp (const FArrayBox& crse,
                                int              crse_comp,
                                FArrayBox&       fine,
                                int              fine_comp,
                                int              ncomp,
                                const Box&       fine_region,
                                const IntVect&   ratio,
                                const Geometry&  crse_geom,
                                const Geometry&  fine_geom,
                                Vector<BCRec> const& bcr,
                                int              /*actual_comp*/,
                                int              /*actual_state*/,
                                RunOn             runon)
{
    BL_PROFILE("CellConservativeLinear::interp()");
    AMREX_ASSERT(bcr.size() >= ncomp);

    AMREX_ASSERT(fine.box().contains(fine_region));

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Box const& cdomain = crse_geom.Domain();
    amrex::ignore_unused(fine_geom);

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();

    const Box& crse_region = CoarseBox(fine_region,ratio);
    Box cslope_bx(crse_region);
    for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
        if (ratio[dim] > 1) {
            cslope_bx.grow(dim,-1);
        }
    }

    FArrayBox ccfab(cslope_bx, ncomp*AMREX_SPACEDIM);
    Array4<Real> const& tmp = ccfab.array();
    Array4<Real const> const& ctmp = ccfab.const_array();

#ifdef AMREX_USE_GPU
    AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
    BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();

    Elixir cceli;
    if (run_on_gpu) { cceli = ccfab.elixir(); }
#else
    BCRec const* bcrp = bcr.data();
#endif

#if (AMREX_SPACEDIM == 1)
    if (crse_geom.IsSPHERICAL()) {
        Real drf = fine_geom.CellSize(0);
        Real rlo = fine_geom.Offset(0);
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                amrex::ignore_unused(j,k);
                mf_cell_cons_lin_interp_mcslope_sph(i, n, tmp, crsearr, crse_comp, ncomp,
                                                    cdomain, ratio, bcrp, drf, rlo);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            amrex::ignore_unused(j,k);
            mf_cell_cons_lin_interp_sph(i, n, finearr, fine_comp, ctmp,
                                        crsearr, crse_comp, ncomp, ratio, drf, rlo);
        });
    } else
#elif (AMREX_SPACEDIM == 2)
    if (crse_geom.IsSPHERICAL()) {
        Real drf = fine_geom.CellSize(0);
        Real dtf = fine_geom.CellSize(1);
        Real rlo = fine_geom.Offset(0);
        Real tlo = fine_geom.Offset(1);
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                amrex::ignore_unused(k);
                mf_cell_cons_lin_interp_mcslope_sph(i, j, n, tmp, crsearr, crse_comp, ncomp,
                                                    cdomain, ratio, bcrp, drf, rlo, dtf, tlo);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            amrex::ignore_unused(k);
            mf_cell_cons_lin_interp_sph(i, j, n, finearr, fine_comp, ctmp, crsearr, crse_comp,
                                        ncomp, ratio, drf, rlo, dtf, tlo);
        });
    } else if (crse_geom.IsRZ()) {
        Real drf = fine_geom.CellSize(0);
        Real rlo = fine_geom.Offset(0);
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                amrex::ignore_unused(k);
                mf_cell_cons_lin_interp_mcslope_rz(i, j, n, tmp, crsearr, crse_comp, ncomp,
                                                   cdomain, ratio, bcrp, drf, rlo);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            amrex::ignore_unused(k);
            mf_cell_cons_lin_interp_rz(i, j, n, finearr, fine_comp, ctmp,
                                       crsearr, crse_comp, ncomp, ratio, drf, rlo);
        });
    } else
#endif
    {
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                mf_cell_cons_lin_interp_mcslope(i,j,k,n, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            mf_cell_cons_lin_interp(i,j,k,n, finearr, fine_comp, ctmp,
                                    crsearr, crse_comp, ncomp, ratio);
        });
    }
}

Box
CellQuadratic::CoarseBox (const Box&     fine,
                          const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellQuadratic::CoarseBox (const Box& fine,
                          int        ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
CellQuadratic::interp (const FArrayBox& crse,
                       int              crse_comp,
                       FArrayBox&       fine,
                       int              fine_comp,
                       int              ncomp,
                       const Box&       fine_region,
                       const IntVect&   ratio,
                       const Geometry&  crse_geom,
                       const Geometry&  fine_geom,
                       Vector<BCRec> const&  bcr,
                       int              /* actual_comp */,
                       int              /* actual_state */,
                       RunOn            runon)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(crse,crse_comp,fine,fine_comp,ncomp,fine_region,
                         ratio,crse_geom,fine_geom,bcr,runon);
    amrex::Abort("1D CellQuadratic::interp not supported");
#else

    BL_PROFILE("CellQuadratic::interp()");
    AMREX_ASSERT(bcr.size() >= ncomp);

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    // Make Box for slopes.
    Box cslope_bx = amrex::coarsen(target_fine_region,ratio);
    AMREX_ASSERT(crse.box().contains(cslope_bx));

    // Are we running on GPU?
    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    // Set up domain for coarse geometry
    Box const& cdomain = crse_geom.Domain();

    // Set up temporary fab (with elixir, as needed) for coarse grid slopes
#if (AMREX_SPACEDIM == 2)
    int nslp = 5; // x, y, x^2, y^2, xy, in that order.
#else  /* AMREX_SPACEDIM == 3 */
    int nslp = 9; // x, y, z, x^2, y^2, z^2, xy, xz, yz, in that order.
#endif /* AMREX_SPACEDIM == 2 */
    FArrayBox sfab(cslope_bx, nslp*ncomp);

#ifdef AMREX_USE_GPU
    // Set up AsyncArray for boundary conditions
    AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
    BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();

    Elixir seli;
    if (run_on_gpu) { seli = sfab.elixir(); }
#else
    BCRec const* bcrp = bcr.data();
#endif

    // Extract pointers to fab data
    Array4<Real>       const&   finearr = fine.array();
    Array4<Real const> const&   crsearr = crse.const_array();
    Array4<Real>       const&  slopearr = sfab.array();
    Array4<Real const> const& cslopearr = sfab.const_array();

    // Compute slopes.
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
    {
        mf_cell_quadratic_calcslope(i, j, k, n,
                                    crsearr, crse_comp,
                                    slopearr,
                                    cdomain, bcrp);
    });

#if (AMREX_SPACEDIM == 2)
    if (crse_geom.IsRZ()) {

        // Get coarse and fine geometry data.
        GeometryData const& cs_geomdata = crse_geom.data();
        GeometryData const& fn_geomdata = fine_geom.data();

        // Compute fine correction.
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                               i, j, k, n,
        {
            mf_cell_quadratic_interp_rz(i, j, k, n,
                                        finearr, fine_comp,
                                        crsearr, crse_comp,
                                        cslopearr,
                                        ratio,
                                        cs_geomdata, fn_geomdata);
        });

    } else { /* crse_geom.IsCartesian() */
#endif /* AMREX_SPACEDIM == 2 */

        // No need for fine geometry data if using Cartesian coordinates.
        amrex::ignore_unused(fine_geom);

        // Compute fine correction.
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                               i, j, k, n,
        {
            mf_cell_quadratic_interp(i, j, k, n,
                                     finearr, fine_comp,
                                     crsearr, crse_comp,
                                     cslopearr,
                                     ratio);
        });

#if (AMREX_SPACEDIM == 2)
    } // geom
#endif /* AMREX_SPACEDIM == 2 */

#endif /*(AMREX_SPACEDIM == 1)*/
}

Box
PCInterp::CoarseBox (const Box& fine,
                     int        ratio)
{
    return amrex::coarsen(fine,ratio);
}

Box
PCInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    return amrex::coarsen(fine,ratio);
}

void
PCInterp::interp (const FArrayBox& crse,
                  int              crse_comp,
                  FArrayBox&       fine,
                  int              fine_comp,
                  int              ncomp,
                  const Box&       fine_region,
                  const IntVect&   ratio,
                  const Geometry& /*crse_geom*/,
                  const Geometry& /*fine_geom*/,
                  Vector<BCRec> const& /*bcr*/,
                  int               /*actual_comp*/,
                  int               /*actual_state*/,
                  RunOn             runon)
{
    BL_PROFILE("PCInterp::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
    {
        amrex::pcinterp_interp(tbx,finearr,fine_comp,ncomp,crsearr,crse_comp,ratio);
    });
}

CellConservativeProtected::CellConservativeProtected ()
    : CellConservativeLinear(true) {}

void
CellConservativeProtected::protect (const FArrayBox& /*crse*/,
                                    int              /*crse_comp*/,
                                    FArrayBox&       fine,
                                    int              /*fine_comp*/,
                                    FArrayBox&       fine_state,
                                    int              /*state_comp*/,
                                    int              ncomp,
                                    const Box&       fine_region,
                                    const IntVect&   ratio,
                                    const Geometry&  crse_geom,
                                    const Geometry&  fine_geom,
                                    Vector<BCRec>&   /*bcr*/,
                                    RunOn            runon)
{
    AMREX_ALWAYS_ASSERT(ratio.allGT(1));

#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(fine,fine_state,
                         ncomp,fine_region,ratio,
                         crse_geom,fine_geom,runon);
    amrex::Abort("1D CellConservativeProtected::protect not supported");
#else
    BL_PROFILE("CellConservativeProtected::protect()");

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    //
    // crse_bx is coarsening of target_fine_region, grown by 1.
    //
    Box crse_bx = CoarseBox(target_fine_region,ratio);

    //
    // cs_bx is coarsening of target_fine_region.
    //
    Box cs_bx(crse_bx);
    cs_bx.grow(-1);

#if (AMREX_SPACEDIM == 2)
    /*
     * Get coarse and fine geometry data.
     */
    GeometryData cs_geomdata = crse_geom.data();
    GeometryData fn_geomdata = fine_geom.data();
#else
    amrex::ignore_unused(crse_geom, fine_geom);
#endif

    // Extract box from fine fab
    const Box& fnbx = fine.box();

    // Extract pointers to fab data
    Array4<Real>       const&   fnarr = fine.array();
    Array4<Real const> const& fnstarr = fine_state.const_array();

    /*
     * Loop over coarse indices.
     */
#if (AMREX_SPACEDIM == 2)
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cs_bx, ic, jc, kc,
    {
        ccprotect_2d(ic, jc, kc, ncomp,
                     fnbx, ratio,
                     cs_geomdata, fn_geomdata,
                     fnarr, fnstarr);
    }); // cs_bx
#else
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cs_bx, ic, jc, kc,
    {
        ccprotect_3d(ic, jc, kc, ncomp,
                     fnbx, ratio,
                     fnarr, fnstarr);
    }); // cs_bx
#endif

#endif /*(AMREX_SPACEDIM == 1)*/

}

Box
CellConservativeQuartic::CoarseBox (const Box& fine,
                                    int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    crse.grow(2);
    return crse;
}

Box
CellConservativeQuartic::CoarseBox (const Box&     fine,
                                    const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

void
CellConservativeQuartic::interp (const FArrayBox&  crse,
                                 int               crse_comp,
                                 FArrayBox&        fine,
                                 int               fine_comp,
                                 int               ncomp,
                                 const Box&        fine_region,
                                 const IntVect&    ratio,
                                 const Geometry&   /* crse_geom */,
                                 const Geometry&   /* fine_geom */,
                                 Vector<BCRec> const& /*bcr*/,
                                 int               /* actual_comp */,
                                 int               /* actual_state */,
                                 RunOn             runon)
{
    BL_PROFILE("CellConservativeQuartic::interp()");
    AMREX_ASSERT(ratio == 2);
    amrex::ignore_unused(ratio);

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    // Extract pointers to fab data
    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp, i, j, k, n,
    {
        ccquartic_interp(i, j, k, n,
                         crsearr, finearr);
    });
}

Box
FaceDivFree::CoarseBox (const Box& fine,
                        int        ratio)
{
    Box crse = amrex::coarsen(fine,ratio).grow(1);
    return crse;
}

Box
FaceDivFree::CoarseBox (const Box&     fine,
                        const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio).grow(1);
    return crse;
}

void
FaceDivFree::interp (const FArrayBox&  /*crse*/,
                     int               /*crse_comp*/,
                     FArrayBox&        /*fine*/,
                     int               /*fine_comp*/,
                     int               /*ncomp*/,
                     const Box&        /*fine_region*/,
                     const IntVect&    /*ratio*/,
                     const Geometry&   /*crse_geom*/,
                     const Geometry&   /*fine_geom*/,
                     Vector<BCRec> const& /*bcr*/,
                     int               /*actual_comp*/,
                     int               /*actual_state*/,
                     RunOn             /*runon*/)
{
    amrex::Abort("FaceDivFree does not work on a single FArrayBox. Call 'interp_arr' instead.");
}

void
FaceDivFree::interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                         const int         crse_comp,
                         Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                         const int         fine_comp,
                         const int         ncomp,
                         const Box&        fine_region,
                         const IntVect&    ratio,
                         Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                         const Geometry&   /*crse_geom */,
                         const Geometry&   fine_geom,
                         Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                         const int         /*actual_comp*/,
                         const int         /*actual_state*/,
                         const RunOn       runon)
{
    BL_PROFILE("FaceDivFree::interp()");

    Array<IndexType, AMREX_SPACEDIM> types;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        { types[d].set(d); }

    // This is currently only designed for octree, where ratio = 2.
    AMREX_ALWAYS_ASSERT(ratio == 2);

    const Box c_fine_region = amrex::coarsen(fine_region, ratio);
    GpuArray<Real, AMREX_SPACEDIM> cell_size = fine_geom.CellSizeArray();

    GpuArray<Array4<const Real>, AMREX_SPACEDIM> crsearr;
    GpuArray<Array4<Real>, AMREX_SPACEDIM> finearr;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> maskarr;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        crsearr[d] = crse[d]->const_array(crse_comp);
        finearr[d] = fine[d]->array(fine_comp);
        if (solve_mask[d] != nullptr)
            { maskarr[d] = solve_mask[d]->const_array(0); }
    }

    // Fuse the launches, 1 for each dimension, into a single launch.
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(c_fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 0,
                                                            crsearr[0], finearr[0], maskarr[0], ratio);
                      }
                  });
              },
              amrex::convert(c_fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 1,
                                                            crsearr[1], finearr[1], maskarr[1], ratio);
                      }
                  });
              },
              amrex::convert(c_fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 2,
                                                            crsearr[2], finearr[2], maskarr[2], ratio);
                      }
                  });
              });

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,c_fine_region,ncomp,i,j,k,n,
    {
        amrex::facediv_int<Real>(i, j, k, fine_comp+n, finearr, ratio, cell_size);
    });
}

Box
CellQuartic::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(3);
    return crse;
}

Box
CellQuartic::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(3);
    return crse;
}

// void
// CellQuartic::interp (const FArrayBox& crse,
//                      int              crse_comp,
//                      FArrayBox&       fine,
//                      int              fine_comp,
//                      int              ncomp,
//                      const Box&       fine_region,
//                      const IntVect&   ratio,
//                      const Geometry&  /*crse_geom*/,
//                      const Geometry&  /*fine_geom*/,
//                      Vector<BCRec> const&  /*bcr*/,
//                      int              /* actual_comp */,
//                      int              /* actual_state */,
//                      RunOn            runon)
// {
//     BL_PROFILE("CellQuartic::interp()");
//     amrex::ignore_unused(ratio);
//     AMREX_ASSERT(ratio == 2);

//     Box target_fine_region = fine_region & fine.box();

//     bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
//     amrex::ignore_unused(run_on_gpu);

//     Array4<Real const> const& crsearr = crse.const_array(crse_comp);
//     Array4<Real>       const& finearr = fine.array(fine_comp);

// #if (AMREX_SPACEDIM == 3)
//     Box bz = amrex::coarsen(target_fine_region, IntVect(2,2,1));
//     bz.grow(IntVect(2,2,0));
//     FArrayBox tmpz(bz, ncomp);
// #ifdef AMREX_USE_GPU
//     Elixir tmpz_eli;
//     if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
// #endif
//     Array4<Real> const& tmpzarr = tmpz.array();
//     AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
//     {
//         cell_quartic_interp_z(i,j,k,n,tmpzarr,crsearr);
//     });
// #endif

// #if (AMREX_SPACEDIM >= 2)
//     Box by = amrex::coarsen(target_fine_region, IntVect(AMREX_D_DECL(2,1,1)));
//     by.grow(IntVect(AMREX_D_DECL(2,0,0)));
//     FArrayBox tmpy(by, ncomp);
// #ifdef AMREX_USE_GPU
//     Elixir tmpy_eli;
//     if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
// #endif
//     Array4<Real> const& tmpyarr = tmpy.array();
// #if (AMREX_SPACEDIM == 2)
//     Array4<Real const> srcarr = crsearr;
// #else
//     Array4<Real const> srcarr = tmpz.const_array();
// #endif
//     AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
//     {
//         cell_quartic_interp_y(i,j,k,n,tmpyarr,srcarr);
//     });
// #endif

// #if (AMREX_SPACEDIM == 1)
//     Array4<Real const> srcarr = crsearr;
// #else
//     srcarr = tmpy.const_array();
// #endif
//     AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
//                                            i, j, k, n,
//     {
//         cell_quartic_interp_x(i,j,k,n,finearr,srcarr);
//    });
// }


void
CellQuartic::interp (const FArrayBox& crse,
                     int              crse_comp,
                     FArrayBox&       fine,
                     int              fine_comp,
                     int              ncomp,
                     const Box&       fine_region,
                     const IntVect&   ratio,
                     const Geometry&  /*crse_geom*/,
                     const Geometry&  /*fine_geom*/,
                     Vector<BCRec> const&  /*bcr*/,
                     int              /* actual_comp */,
                     int              /* actual_state */,
                     RunOn            runon)
{
    BL_PROFILE("5th-order Lagrange::interp()");
    // amrex::Print() << "5th-order Lagrange interp: fine_box " << fine.box() << "\n";
    // amrex::Print() << "5th-order Lagrange interp: fine_region " << fine_region << "\n";
    amrex::ignore_unused(ratio);
    AMREX_ASSERT(ratio == 2);

    Box target_fine_region = fine_region & fine.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::coarsen(target_fine_region, IntVect(2,2,1));
    bz.grow(IntVect(3,3,0));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        Lagrange5_interp_z(i,j,k,n,tmpzarr,crsearr);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::coarsen(target_fine_region, IntVect(AMREX_D_DECL(2,1,1)));
    by.grow(IntVect(AMREX_D_DECL(3,0,0)));
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = crsearr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif
    // amrex::Print() << "crse.box()" << crse.box() << std::endl;

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
    {
        Lagrange5_interp_y(i,j,k,n,tmpyarr,srcarr);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = crsearr;
#else
    srcarr = tmpy.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                           i, j, k, n,
    {
        Lagrange5_interp_x(i,j,k,n,finearr,srcarr);
   });
}

void
CellQuartic::restrict (const FArrayBox& fine,
                     int              fine_comp,
                     FArrayBox&       crse,
                     int              crse_comp,
                     int              ncomp,
                     const Box&       crse_region,
                     const IntVect&   ratio,
                     const Geometry&  /*fine_geom*/,
                     const Geometry&  /*crse_geom*/,
                     Vector<BCRec> const&  /*bcr*/,
                     int              /* actual_comp */,
                     int              /* actual_state */,
                     RunOn            runon)
{
    BL_PROFILE("CellQuartic::restrict()");
    amrex::ignore_unused(ratio);
    AMREX_ASSERT(ratio == 2);

    Box target_crse_region = crse_region & crse.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& finearr = fine.const_array(fine_comp);
    Array4<Real>       const& crsearr = crse.array(crse_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::refine(target_crse_region, IntVect(2,2,1));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        cell_quartic_restrict_z(i,j,k,n,tmpzarr,finearr);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::refine(target_crse_region, IntVect(AMREX_D_DECL(2,1,1)));
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = finearr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
    {
        cell_quartic_restrict_y(i,j,k,n,tmpyarr,srcarr);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = finearr;
#else
    srcarr = tmpy.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, ncomp,
                                           i, j, k, n,
    {
        cell_quartic_restrict_x(i,j,k,n,crsearr,srcarr);
    });
}


Box
CellWENO::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(3);
    return crse;
}

Box
CellWENO::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(3);
    return crse;
}

void
CellWENO::interp (const FArrayBox& crse,
                     int              crse_comp,
                     FArrayBox&       fine,
                     int              fine_comp,
                     int              ncomp,
                     const Box&       fine_region,
                     const IntVect&   ratio,
                     const Geometry&  /*crse_geom*/,
                     const Geometry&  /*fine_geom*/,
                     Vector<BCRec> const&  /*bcr*/,
                     int              /* actual_comp */,
                     int              /* actual_state */,
                     RunOn            runon)
{
    BL_PROFILE("CellWENO::interp()");

    // Support both ratio=2 and ratio=4
    AMREX_ASSERT(ratio == 2 || ratio == 4);

    Box target_fine_region = fine_region & fine.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::coarsen(target_fine_region, IntVect(ratio[0],ratio[1],1));
    bz.grow(IntVect(3,3,0));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        central_weno_interp_z(i,j,k,n,tmpzarr,crsearr,ratio);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::coarsen(target_fine_region, IntVect(AMREX_D_DECL(ratio[0],1,1)));
    by.grow(IntVect(AMREX_D_DECL(3,0,0)));
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = crsearr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
    {
        central_weno_interp_y(i,j,k,n,tmpyarr,srcarr,ratio);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = crsearr;
#else
    srcarr = tmpy.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                           i, j, k, n,
    {
        central_weno_interp_x(i,j,k,n,finearr,srcarr,ratio);
   });
}


void
CellWENO::restrict (const FArrayBox& fine,
                    int              fine_comp,
                    FArrayBox&       crse,
                    int              crse_comp,
                    int              ncomp,
                    const Box&       crse_region,
                    const IntVect&   ratio,
                    const Geometry&  /*fine_geom*/,
                    const Geometry&  /*crse_geom*/,
                    Vector<BCRec> const&  /*bcr*/,
                    int              /* actual_comp */,
                    int              /* actual_state */,
                    RunOn            runon)
{
    BL_PROFILE("CellWENO::restrict()");
    
    // Support both ratio=2 and ratio=4
    AMREX_ASSERT(ratio == 2 || ratio == 4);

    Box target_crse_region = crse_region & crse.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& finearr = fine.const_array(fine_comp);
    Array4<Real>       const& crsearr = crse.array(crse_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::refine(target_crse_region, IntVect(ratio[0],ratio[1],1));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        weno_restrict_z(i,j,k,n,tmpzarr,finearr,ratio);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::refine(target_crse_region, IntVect(AMREX_D_DECL(ratio[0],1,1)));
    by.grow(IntVect(AMREX_D_DECL(ratio[0],0,0))); // halo for x-stage i±1 stencil
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = finearr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
    {
        weno_restrict_y(i,j,k,n,tmpyarr,srcarr,ratio);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = finearr;
#else
    srcarr = tmpy.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, ncomp,
                                           i, j, k, n,
    {
        weno_restrict_x(i,j,k,n,crsearr,srcarr,ratio);
    });
}



Box
Mortar2D::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    return crse;
}

Box
Mortar2D::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    return crse;
}

void
Mortar2D::interp (const FArrayBox& crse,
                    int              crse_comp,
                    FArrayBox&       fine,
                    int              fine_comp,
                    int              ncomp,
                    const Box&       fine_region,
                    const IntVect&   ratio,
                    const Geometry&  /*crse_geom*/,
                    const Geometry&  /*fine_geom*/,
                    Vector<BCRec> const&  /*bcr*/,
                    int              /* actual_comp */,
                    int              /* actual_state */,
                    RunOn            runon)
{
    BL_PROFILE("Mortar2D::interp()");
    
    // Support both ratio=2 and ratio=4
    AMREX_ASSERT(ratio == 2 || ratio == 4);

    Box target_fine_region = fine_region & fine.box();

    // amrex::Print() << "target_fine_region=" << target_fine_region   << std::endl;

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    // amrex::Print() << "crse_comp= " << crse_comp
    //                << "\nfine_comp= " << fine_comp << std::endl;
    // amrex::Abort("aa");

    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 2)
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp/sd_space, i, j, k, n,
    {
        mortar_interp(i,j,k,n,finearr,crsearr,ratio);
    });
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Mortar2D::mortar_interp(const int i, const int j, const int k, const int n, 
                Array4<Real> const& finearr, Array4<const Real> const& crsearr,
                const IntVect&   ratio)
{
    if (Mortar2D::type == Type::OrderRef)
    {
        // define a smaller and a larger 2D matrix 
        Real uc[sd_order][sd_order];
        Real uf[sd_order*2][sd_order*2]; // row <i> of u corresponds to physical axis <x>


        int ic = amrex::coarsen(i, ratio[0]);
        int jc = amrex::coarsen(j, ratio[1]);
        int kc = (AMREX_SPACEDIM>2)? amrex::coarsen(k, ratio[2]) : 0;

        int nc;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // 5×5 coarse state
                nc = box2point(s, m, 0, n);
                uc[s][m] = crsearr(ic,jc,kc,nc);
            }
        }

        Real tmp[sd_order][sd_order*2];
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order*2; ++m) {
                tmp[s][m] = Real(0.0);
                for (int q = 0; q < sd_order; ++q) {
                    tmp[s][m] += uc[s][q] * Py_prol[q][m];
                }
            }
        }

        for (int s = 0; s < sd_order*2; ++s) {
            for (int m = 0; m < sd_order*2; ++m) {
                uf[s][m] = Real(0.0);
                for (int q = 0; q < sd_order; ++q) {
                    uf[s][m] += Px_prol[s][q] * tmp[q][m];
                }
            }
        }

        int ioff = i-ratio[0]*ic, joff = j-ratio[1]*jc;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                nc = box2point(s, m, 0, n);
                finearr(i,j,k,nc) = uf[ioff*sd_order+s][joff*sd_order+m];
            }
        }
        
        return;
    } else if (Mortar2D::type == Type::ScaleRef)
    {
        // define four 2D mass matrices with equivalent sizes
        Real uc[sd_order][sd_order];

        int ic = amrex::coarsen(i, ratio[0]);
        int jc = amrex::coarsen(j, ratio[1]);
        int kc = (AMREX_SPACEDIM>2)? amrex::coarsen(k, ratio[2]) : 0;

        // amrex::Print() << "ic=" << ic << " jc=" << jc << std::endl;

        int nc;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // 5*5 crse state
                nc = box2point(s, m, 0, n);
                uc[s][m] = crsearr(ic,jc,kc,nc);
            }
        }

        Real tmp[sd_order][sd_order];
    
        int joff = j-ratio[1]*jc;
        // 5*5 temporary state
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                tmp[s][m] = Real(0);
                for (int q = 0; q < sd_order; ++q) {
                    tmp[s][m] += uc[s][q] * (Real(1.0) * Py2D_prol[joff][q][m]);
                }
            }
        }

        int ioff = i-ratio[0]*ic;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                nc = box2point(s, m, 0, n);
                finearr(i,j,k,nc) = Real(0);
                for (int q = 0; q < sd_order; ++q) {
                    finearr(i,j,k,nc) += (Real(1.0) * Px2D_prol[ioff][s][q]) * tmp[q][m];
                }
            }
        }
        // ---------------- debug ----------------------------------
        // for (int z = 0; z < sd_order; ++z) {
        //     amrex::Print() << "\nic= " << ic << "jc= " << jc << "\ncoarse fab=" << crsearr(ic,jc,kc,z) << " ";
        // }
        // for (int z = 0; z < sd_order; ++z) {
        //     amrex::Print() << "\ni= " << i << "j= " << j << "\nfine fab=" << finearr(i,j,k,z) << " ";
        // }
        return;

    } else {
        amrex::Abort("Unknown strategy for c->f prolongation");
    }
}


void
Mortar2D::restrict (const FArrayBox& fine,
                    int              fine_comp,
                    FArrayBox&       crse,
                    int              crse_comp,
                    int              ncomp,
                    const Box&       crse_region,
                    const IntVect&   ratio,
                    const Geometry&  /*fine_geom*/,
                    const Geometry&  /*crse_geom*/,
                    Vector<BCRec> const&  /*bcr*/,
                    int              /* actual_comp */,
                    int              /* actual_state */,
                    RunOn            runon)
{
    BL_PROFILE("Mortar2D::restrict()");
    
    // Support both ratio=2 and ratio=4
    AMREX_ASSERT(ratio == 2 || ratio == 4);

    Box target_crse_region = crse_region & crse.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& finearr = fine.const_array(fine_comp);
    Array4<Real>       const& crsearr = crse.array(crse_comp);


#if (AMREX_SPACEDIM == 2)
    auto const& destarr = crsearr;
    auto const& srcarr  = finearr;
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, ncomp/sd_space, i, j, k, n,
    {
        mortar_restrict(i,j,k,n,destarr,srcarr,ratio);
        // mortar_restrict_flatten(i,j,k,n,destarr,srcarr,ratio);
    });
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Mortar2D::mortar_restrict(const int i, const int j, const int k, const int n, 
                Array4<Real> const crsearr, Array4<const Real> const& finearr,
                const IntVect&   ratio)
{
    AMREX_ASSERT(ratio[0]==2 && ratio[1]==2);
    if (Mortar2D::type == Type::OrderRef)
    {
        // define a smaller and a larger 2D matrix 
        Real uc[sd_order][sd_order] = {};
        Real uf[sd_order*2][sd_order*2] = {};


        int ii = i*ratio[0];
        int jj = j*ratio[1];
        int kk = (AMREX_SPACEDIM>2)? k*ratio[2] : 0;

        // assemble larger 2D working matrix on finer grids
        int nc;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // 10×10 fine state
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    for (int joff = 0; joff < ratio[1]; ++joff) {
                        nc = box2point(s, m, 0, n);
                        uf[ioff*sd_order+s][joff*sd_order+m] = finearr(ii+ioff,jj+joff,kk,nc);
                    }
                }
            }
        }

        // execuate y-direction projection first
        Real tmp[sd_order*2][sd_order] = {};
        for (int s = 0; s < sd_order*2; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // tmp[s][m] = Real(0.0);
                for (int q = 0; q < sd_order*2; ++q) {
                    tmp[s][m] += uf[s][q] * Py[q][m];
                }
            }
        }

        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // uc[s][m] = Real(0.0);
                for (int q = 0; q < sd_order*2; ++q) {
                    uc[s][m] += Px[s][q]*tmp[q][m];
                }
                // 5×5 coarse state
                nc = box2point(s, m, 0, n);
                crsearr(i,j,k,nc) = uc[s][m];
            }
        }
        return;

    } else if (Mortar2D::type == Type::ScaleRef)
    {
        // define four 2D mass matrices with equivalent sizes
        Real uc[sd_order][sd_order];
        Real uf[2][2][sd_order][sd_order]; // row <i> of u corresponds to physical axis <x>


        int ii = i*ratio[0];
        int jj = j*ratio[1];
        int kk = (AMREX_SPACEDIM>2)? k*ratio[2] : 0;

        for (int ioff = 0; ioff < ratio[0]; ++ioff) {
            for (int joff = 0; joff < ratio[1]; ++joff) {
                // 5*5 fine state
                for (int s = 0; s < sd_order; ++s) {
                    for (int m = 0; m < sd_order; ++m) {
                        int nc = box2point(s, m, 0, n);
                        uf[ioff][joff][s][m] = finearr(ii+ioff,jj+joff,kk,nc);
                    }
                }
            }
        }

        Real tmp[2][2][sd_order][sd_order];
        // 5*5 fine state
        for (int ioff = 0; ioff < ratio[0]; ++ioff) {
            for (int joff = 0; joff < ratio[1]; ++joff) {
                for (int s = 0; s < sd_order; ++s) {
                    for (int m = 0; m < sd_order; ++m) {
                        tmp[ioff][joff][s][m] = Real(0);
                        for (int q = 0; q < sd_order; ++q) {
                            tmp[ioff][joff][s][m] += uf[ioff][joff][s][q] * Py2D[joff][q][m];
                        }
                    }
                }
            }
        }

        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                uc[s][m] = Real(0);
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    for (int joff = 0; joff < ratio[1]; ++joff) {
                        for (int q = 0; q < sd_order; ++q) {
                            uc[s][m] += Px2D[ioff][s][q] * tmp[ioff][joff][q][m];
                        }
                    }
                }
                int nc = box2point(s, m, 0, n);
                crsearr(i,j,k,nc) = uc[s][m];
            }
        }
        return;

    } else {
        amrex::Abort("Unknown strategy for f->c restriction");
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Mortar2D::mortar_restrict_flatten(const int i, const int j, const int k, const int n, 
                                  Array4<Real> const crsearr, Array4<const Real> const& finearr,
                                  const IntVect&   ratio)
{
    AMREX_ASSERT(ratio[0]==2 && ratio[1]==2);
    if (Mortar2D::type == Type::OrderRef)
    {
        // define the coarse- and fine-state vectors 
        int fsd_order = ratio[0]*sd_order;
        int fsd_space = fsd_order*fsd_order;
        Real uc[sd_space]  = {};
        Real uf[fsd_space] = {};

        // indices on the fine grids
        int ii = i*ratio[0];
        int jj = j*ratio[1];
        int kk = (AMREX_SPACEDIM>2)? k*ratio[2] : 0;

        // 1D flattening
        int nc, ixf, iyf;
        for (int iy = 0; iy < sd_order; ++iy) {
            for (int ix = 0; ix < sd_order; ++ix) {
                // within a big coarse cell
                for (int joff = 0; joff < ratio[1]; ++joff) {
                    for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                        nc = box2point(ix, iy, 0, n);
                        ixf = ix + ioff*sd_order; // begin from 0
                        iyf = iy + joff*sd_order; // begin from 0
                        int kf = fsd_order*ixf + iyf;
                        uf[kf] = finearr(ii+ioff,jj+joff,kk,nc);
                    }
                }
            }
        }

        // fine-to-coarse rojection
        for (int ix = 0; ix < sd_space; ++ix) {
            for (int iy = 0; iy < fsd_space; ++iy) {
                uc[ix] += P2D[ix][iy]*uf[iy];
            }
        }

        // Unflattening
        for (int iy = 0; iy < sd_order; ++iy) {
            for (int ix = 0; ix < sd_order; ++ix) {
                int kc = sd_order*ix + iy;
                nc = box2point(ix, iy, 0, n);
                crsearr(i,j,k,nc) = uc[kc];
            }
        }
        return;

    } else if (Mortar2D::type == Type::ScaleRef)
    {
        // define four 2D mass matrices with equivalent sizes
        Real uc[sd_order][sd_order];
        Real uf[2][2][sd_order][sd_order]; // row <i> of u corresponds to physical axis <x>


        int ii = i*ratio[0];
        int jj = j*ratio[1];
        int kk = (AMREX_SPACEDIM>2)? k*ratio[2] : 0;

        int nc;
        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                // 2D 5*5 fine state
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    for (int joff = 0; joff < ratio[1]; ++joff) {
                        nc = box2point(s, m, 0, n);
                        uf[ioff][joff][s][m] = finearr(ii+ioff,jj+joff,kk,nc);
                    }
                }
            }
        }

        Real tmp[2][2][sd_order][sd_order];
        // 5*5 fine state
        for (int ioff = 0; ioff < ratio[0]; ++ioff) {
            for (int joff = 0; joff < ratio[1]; ++joff) {
                for (int s = 0; s < sd_order; ++s) {
                    for (int m = 0; m < sd_order; ++m) {
                        tmp[ioff][joff][s][m] = Real(0);
                        for (int q = 0; q < sd_order; ++q) {
                            tmp[ioff][joff][s][m] += uf[ioff][joff][s][q] * Py2D[joff][q][m];
                        }
                    }
                }
            }
        }

        for (int s = 0; s < sd_order; ++s) {
            for (int m = 0; m < sd_order; ++m) {
                uc[s][m] = Real(0);
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    for (int joff = 0; joff < ratio[1]; ++joff) {
                        for (int q = 0; q < sd_order; ++q) {
                            uc[s][m] += Px2D[ioff][s][q] * tmp[ioff][joff][q][m];
                        }
                    }
                }
                nc = box2point(s, m, 0, n);
                crsearr(i,j,k,nc) = uc[s][m];
            }
        }
        return;

    } else {
        amrex::Abort("Unknown strategy for f->c restriction");
    }
}

Box
HermiteWENO2D::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(IntVect(AMREX_D_DECL(1,1,0)));
    return crse;
}

Box
HermiteWENO2D::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(IntVect(AMREX_D_DECL(1,1,0)));
    return crse;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HermiteWENO2D::hweno_interp_z (int i, int j, int k, int n,
                               Array4<Real> const& tmparr,
                               Array4<Real const> const& srcarr,
                               IntVect const& ratio,
                               Real hdir) noexcept
{
    // In d-by-d interp(), z-stage loops over `bz`, where x- and y-indices are already coarse.
    const int ic = i;
    const int jc = j;
    const int kc = amrex::coarsen(k, ratio[2]);

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HermiteWENO2D::hweno_interp_y (int i, int j, int k, int n,
                               Array4<Real> const& tmparr,
                               Array4<Real const> const& srcarr,
                               IntVect const& ratio,
                               Real hdir,
                               Real eos_gamma) noexcept
{
    if (n != 0) { return; }
    const int nvar = srcarr.nComp() / sd_space_hweno;
    if (nvar <= 0) { return; }

    // In d-by-d interp(), y-stage loops over `by`, where x-index is already coarse.
    const int ic = i;
    const int jc = amrex::coarsen(j, ratio[1]);
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int child = j - ratio[1]*jc;
    const IntVect ivm(AMREX_D_DECL(ic, jc-1, kc));
    const IntVect ivc(AMREX_D_DECL(ic, jc,   kc));
    const IntVect ivp(AMREX_D_DECL(ic, jc+1, kc));

    const SmoothRegion region = (child == 0) ? SmoothRegion::ChildLeft
                                             : SmoothRegion::ChildRight;
    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(line, s, 0, n); // y-dir: line is x-node
            Um[s] = srcarr(ivm[0], ivm[1], ivm[2], off);
            U0[s] = srcarr(ivc[0], ivc[1], ivc[2], off);
            Up[s] = srcarr(ivp[0], ivp[1], ivp[2], off);
        }

        const auto mm = NodalToMoments1D(Um.data(), hdir);
        const auto m0 = NodalToMoments1D(U0.data(), hdir);
        const auto mp = NodalToMoments1D(Up.data(), hdir);
        const auto cand = Build4Candidates(mm.ubar, m0.ubar, mp.ubar,
                                           mm.g1, mp.g1, mm.g2, mp.g2, hdir);

        GpuArray<Real,4> beta{};
        for (int kk = 0; kk < 4; ++kk) {
            beta[kk] = BetaFromCubic(cand[kk], SmoothRegion::Parent);
        }
        const auto omega = ZWeights(beta);

        GpuArray<GpuArray<Real,3>,5> bH_cons{};
        for (int h = 0; h < 5; ++h) {
            for (int m = 0; m < 3; ++m) {
                for (int a = 0; a < 5; ++a) {
                    bH_cons[h][m] += R[h][a] * bH_hydro[a][m];
                }
            }
        }

        for (int nv = 0; nv < nvar; ++nv) {
            GpuArray<Real,3> bH{};
            if (hweno_is_hydro_comp(nv)) {
                const int h = (nv == MRHO) ? 0
                            : (nv == MU)   ? 1
                            : (nv == MV)   ? 2
                            : (nv == MW)   ? 3
                                           : 4;
                bH = bH_cons[h];
            } else {
                GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(line, s, 0, nv);
                    Um[s] = srcarr(ivm[0], ivm[1], ivm[2], off);
                    U0[s] = srcarr(ivc[0], ivc[1], ivc[2], off);
                    Up[s] = srcarr(ivp[0], ivp[1], ivp[2], off);
                }
                const auto mm = NodalToMoments1D(Um.data(), hdir);
                const auto m0 = NodalToMoments1D(U0.data(), hdir);
                const auto mp = NodalToMoments1D(Up.data(), hdir);
                const auto cand = Build3Candidates(mm.ubar, m0.ubar, mp.ubar,
                                                   mm.g1, mp.g1, hdir);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], region);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH[m] += omega[kk] * cand[kk][m];
                    }
                }
            }

            GpuArray<Real,sd_order_hweno> up{};
            for (int q = 0; q < sd_order_hweno; ++q) {
                up[q] = EvalQuadratic(bH, xi_sol[q]/2.0);
            }
            for (int s = 0; s < sd_order_hweno; ++s) {
                Real val = Real(0.0);
                for (int q = 0; q < sd_order_hweno; ++q) {
                    val += up[q] * P1DProjY(child, q, s);
                }
                const int off = box2point(line, s, 0, nv);
                tmparr(i, j, k, off) = val;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HermiteWENO2D::hweno_interp_x (int i, int j, int k, int n,
                               Array4<Real> const& finearr,
                               Array4<Real const> const& srcarr,
                               IntVect const& ratio,
                               Real hdir,
                               Real eos_gamma) noexcept
{
    if (n != 0) { return; }
    const int nvar = srcarr.nComp() / sd_space_hweno;
    if (nvar <= 0) { return; }

    const int ic = amrex::coarsen(i, ratio[0]);
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int child = i - ratio[0]*ic;
    const IntVect ivm(AMREX_D_DECL(ic-1, jc, kc));
    const IntVect ivc(AMREX_D_DECL(ic,   jc, kc));
    const IntVect ivp(AMREX_D_DECL(ic+1, jc, kc));

    const SmoothRegion region = (child == 0) ? SmoothRegion::ChildLeft
                                             : SmoothRegion::ChildRight;
    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(s, line, 0, n); // x-dir: line is y-node
            Um[s] = srcarr(ivm[0], ivm[1], ivm[2], off);
            U0[s] = srcarr(ivc[0], ivc[1], ivc[2], off);
            Up[s] = srcarr(ivp[0], ivp[1], ivp[2], off);
        }

        const auto mm = NodalToMoments1D(Um.data(), hdir);
        const auto m0 = NodalToMoments1D(U0.data(), hdir);
        const auto mp = NodalToMoments1D(Up.data(), hdir);
        const auto cand = Build4Candidates(mm.ubar, m0.ubar, mp.ubar,
                                           mm.g1, mp.g1, mm.g2, mp.g2, hdir);

        GpuArray<Real,4> beta{};
        for (int kk = 0; kk < 4; ++kk) {
            beta[kk] = BetaFromCubic(cand[kk], SmoothRegion::Parent);
        }
        const auto omega = ZWeights(beta);

        GpuArray<GpuArray<Real,3>,5> bH_cons{};
        for (int h = 0; h < 5; ++h) {
            for (int m = 0; m < 3; ++m) {
                for (int a = 0; a < 5; ++a) {
                    bH_cons[h][m] += R[h][a] * bH_hydro[a][m];
                }
            }
        }

        for (int nv = 0; nv < nvar; ++nv) {
            GpuArray<Real,3> bH{};
            if (hweno_is_hydro_comp(nv)) {
                const int h = (nv == MRHO) ? 0
                            : (nv == MU)   ? 1
                            : (nv == MV)   ? 2
                            : (nv == MW)   ? 3
                                           : 4;
                bH = bH_cons[h];
            } else {
                GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(s, line, 0, nv);
                    Um[s] = srcarr(ivm[0], ivm[1], ivm[2], off);
                    U0[s] = srcarr(ivc[0], ivc[1], ivc[2], off);
                    Up[s] = srcarr(ivp[0], ivp[1], ivp[2], off);
                }
                const auto mm = NodalToMoments1D(Um.data(), hdir);
                const auto m0 = NodalToMoments1D(U0.data(), hdir);
                const auto mp = NodalToMoments1D(Up.data(), hdir);
                const auto cand = Build3Candidates(mm.ubar, m0.ubar, mp.ubar,
                                                   mm.g1, mp.g1, hdir);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], region);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH[m] += omega[kk] * cand[kk][m];
                    }
                }
            }

            GpuArray<Real,sd_order_hweno> up{};
            for (int q = 0; q < sd_order_hweno; ++q) {
                up[q] = EvalQuadratic(bH, xi_sol[q]/2.0);
            }
            for (int s = 0; s < sd_order_hweno; ++s) {
                Real val = Real(0.0);
                for (int q = 0; q < sd_order_hweno; ++q) {
                    val += P1DProjX(child, s, q) * up[q];
                }

                const int off = box2point(s, line, 0, nv);
                finearr(i, j, k, off) = val;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HermiteWENO2D::hweno_restrict_y (int i, int j, int k, int n,
                                 Array4<Real> const& tmparr,
                                 Array4<Real const> const& srcarr,
                                 IntVect const& ratio,
                                 Real hfine,
                                 Real hcrse,
                                 Real eos_gamma,
                                 int j_valid_lo,
                                 int j_valid_hi) noexcept
{
    if (n != 0) { return; }
    const int nvar = srcarr.nComp() / sd_space_hweno;
    if (nvar <= 0) { return; }

    const int ic = i;
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;
    const bool has_left  = (jc > j_valid_lo);
    const bool has_right = (jc < j_valid_hi);
    const bool use_central_hweno = (has_left && has_right);
    const bool use_left_biased_p3 = (!has_left && has_right);
    const bool use_right_biased_p2 = (has_left && !has_right);

    const int j0f0 = ratio[1] * jc;
    const int j0f1 = j0f0 + 1;
    const IntVect iv0f0(AMREX_D_DECL(ic, j0f0, kc));
    const IntVect iv0f1(AMREX_D_DECL(ic, j0f1, kc));

    IntVect ivmf0(AMREX_D_DECL(0, 0, 0));
    IntVect ivmf1(AMREX_D_DECL(0, 0, 0));
    IntVect ivpf0(AMREX_D_DECL(0, 0, 0));
    IntVect ivpf1(AMREX_D_DECL(0, 0, 0));
    if (has_left) {
        const int jmf0 = ratio[1] * (jc - 1);
        const int jmf1 = jmf0 + 1;
        ivmf0 = IntVect(AMREX_D_DECL(ic, jmf0, kc));
        ivmf1 = IntVect(AMREX_D_DECL(ic, jmf1, kc));
    }
    if (has_right) {
        const int jpf0 = ratio[1] * (jc + 1);
        const int jpf1 = jpf0 + 1;
        ivpf0 = IntVect(AMREX_D_DECL(ic, jpf0, kc));
        ivpf1 = IntVect(AMREX_D_DECL(ic, jpf1, kc));
    }

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<GpuArray<Real,sd_order_hweno>,5> Um0_h{}, Um1_h{}, U00_h{}, U01_h{}, Up0_h{}, Up1_h{};
        for (int h = 0; h < 5; ++h) {
            const int comp = hweno_hydro_comp(h);
            for (int s = 0; s < sd_order_hweno; ++s) {
                const int off = box2point(line, s, 0, comp);
                U00_h[h][s] = srcarr(iv0f0[0], iv0f0[1], iv0f0[2], off);
                U01_h[h][s] = srcarr(iv0f1[0], iv0f1[1], iv0f1[2], off);
                if (has_left) {
                    Um0_h[h][s] = srcarr(ivmf0[0], ivmf0[1], ivmf0[2], off);
                    Um1_h[h][s] = srcarr(ivmf1[0], ivmf1[1], ivmf1[2], off);
                }
                if (has_right) {
                    Up0_h[h][s] = srcarr(ivpf0[0], ivpf0[1], ivpf0[2], off);
                    Up1_h[h][s] = srcarr(ivpf1[0], ivpf1[1], ivpf1[2], off);
                }
            }
        }

        GpuArray<RestrictInput1D,5> rin_h{};
        if (use_central_hweno) {
            for (int h = 0; h < 5; ++h) {
                const auto m_m5q4 = NodalToMoments1D(Um0_h[h].data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1_h[h].data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00_h[h].data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01_h[h].data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0_h[h].data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1_h[h].data(), hfine);
                rin_h[h] = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                                m_p1q4, m_p3q4, m_p5q4);
            }
        }

        GpuArray<GpuArray<Real,3>,5> bH_cons{};
        bool hydro_ready = false;
        if (use_central_hweno) {
            GpuArray<Real,5> Uj{};
            Uj[0] = rin_h[0].Uj;
            Uj[1] = rin_h[1].Uj;
            Uj[2] = rin_h[2].Uj;
            Uj[3] = rin_h[3].Uj;
            Uj[4] = rin_h[4].Uj;
            GpuArray<GpuArray<Real,5>,5> L{}, R{};
            GpuArray<Real,5> lambda{};
            hweno_build_eigensystem_ideal5(L, R, lambda, Uj, 1, eos_gamma);

            GpuArray<GpuArray<Real,3>,5> bH_hydro{};
            for (int a = 0; a < 5; ++a) {
                Real ujm1 = Real(0.0), uj = Real(0.0), ujp1 = Real(0.0);
                Real gjm1 = Real(0.0), gjp1 = Real(0.0);
                for (int h = 0; h < 5; ++h) {
                    ujm1 += L[a][h] * rin_h[h].Ujm1;
                    uj   += L[a][h] * rin_h[h].Uj;
                    ujp1 += L[a][h] * rin_h[h].Ujp1;
                    gjm1 += L[a][h] * rin_h[h].Gjm1_1;
                    gjp1 += L[a][h] * rin_h[h].Gjp1_1;
                }
                const auto cand = Build3Candidates(ujm1, uj, ujp1, gjm1, gjp1, hcrse);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], SmoothRegion::Parent);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH_hydro[a][m] += omega[kk] * cand[kk][m];
                    }
                }
            }
            for (int h = 0; h < 5; ++h) {
                for (int m = 0; m < 3; ++m) {
                    for (int a = 0; a < 5; ++a) {
                        bH_cons[h][m] += R[h][a] * bH_hydro[a][m];
                    }
                }
            }
            hydro_ready = true;
        }

        for (int nv = 0; nv < nvar; ++nv) {
            GpuArray<Real,3> bH{};
            if (hydro_ready && hweno_is_hydro_comp(nv)) {
                const int h = (nv == MRHO) ? 0
                            : (nv == MU)   ? 1
                            : (nv == MV)   ? 2
                            : (nv == MW)   ? 3
                                           : 4;
                bH = bH_cons[h];
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(line, s, 0, nv);
                    tmparr(i, j, k, off) = EvalQuadratic(bH, xi_sol[s]/2.0);
                }
                continue;
            }

            GpuArray<Real,sd_order_hweno> Um0{}, Um1{}, U00{}, U01{}, Up0{}, Up1{};
            for (int s = 0; s < sd_order_hweno; ++s) {
                const int off = box2point(line, s, 0, nv);
                U00[s] = srcarr(iv0f0[0], iv0f0[1], iv0f0[2], off);
                U01[s] = srcarr(iv0f1[0], iv0f1[1], iv0f1[2], off);
                if (has_left) {
                    Um0[s] = srcarr(ivmf0[0], ivmf0[1], ivmf0[2], off);
                    Um1[s] = srcarr(ivmf1[0], ivmf1[1], ivmf1[2], off);
                }
                if (has_right) {
                    Up0[s] = srcarr(ivpf0[0], ivpf0[1], ivpf0[2], off);
                    Up1[s] = srcarr(ivpf1[0], ivpf1[1], ivpf1[2], off);
                }
            }

            if (use_central_hweno) {
                const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
                const auto rin = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                                      m_p1q4, m_p3q4, m_p5q4);
                const auto cand = Build3Candidates(rin.Ujm1, rin.Uj, rin.Ujp1,
                                                   rin.Gjm1_1, rin.Gjp1_1, hcrse);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], SmoothRegion::Parent);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH[m] += omega[kk] * cand[kk][m];
                    }
                }
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(line, s, 0, nv);
                    tmparr(i, j, k, off) = EvalQuadratic(bH, xi_sol[s]/2.0);
                }
                continue;
            }

            if (use_left_biased_p3) {
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
                const Real Uj = Real(0.5) * (m_m1q4.ubar + m_p1q4.ubar);
                const Real Ujp1 = Real(0.5) * (m_p3q4.ubar + m_p5q4.ubar);
                const Real a = Ujp1 - Uj;
                const Real b = Uj;
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(line, s, 0, nv);
                    const Real xi = xi_sol[s] / Real(2.0);
                    tmparr(i, j, k, off) = a * xi + b;
                }
            } else if (use_right_biased_p2) {
                const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const Real Ujm1 = Real(0.5) * (m_m5q4.ubar + m_m3q4.ubar);
                const Real Uj = Real(0.5) * (m_m1q4.ubar + m_p1q4.ubar);
                const Real a = Uj - Ujm1;
                const Real b = Uj;
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(line, s, 0, nv);
                    const Real xi = xi_sol[s] / Real(2.0);
                    tmparr(i, j, k, off) = a * xi + b;
                }
            } else {
                for (int s = 0; s < sd_order_hweno; ++s) {
                    Real val = Real(0.0);
                    for (int q = 0; q < sd_order_hweno; ++q) {
                        val += P1DRestrictY(0, q, s) * U00[q];
                        val += P1DRestrictY(1, q, s) * U01[q];
                    }
                    const int off = box2point(line, s, 0, nv);
                    tmparr(i, j, k, off) = val;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HermiteWENO2D::hweno_restrict_x (int i, int j, int k, int n,
                                 Array4<Real> const& crsearr,
                                 Array4<Real const> const& srcarr,
                                 IntVect const& ratio,
                                 Real hfine,
                                 Real hcrse,
                                 Real eos_gamma,
                                 int i_valid_lo,
                                 int i_valid_hi) noexcept
{
    if (n != 0) { return; }
    const int nvar = srcarr.nComp() / sd_space_hweno;
    if (nvar <= 0) { return; }

    const int ic = i;
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;
    const bool has_left  = (ic > i_valid_lo);
    const bool has_right = (ic < i_valid_hi);
    const bool use_central_hweno = (has_left && has_right);
    const bool use_left_biased_p3 = (!has_left && has_right);
    const bool use_right_biased_p2 = (has_left && !has_right);

    const int i0f0 = ratio[0] * ic;
    const int i0f1 = i0f0 + 1;
    const IntVect iv0f0(AMREX_D_DECL(i0f0, jc, kc));
    const IntVect iv0f1(AMREX_D_DECL(i0f1, jc, kc));

    IntVect ivmf0(AMREX_D_DECL(0, 0, 0));
    IntVect ivmf1(AMREX_D_DECL(0, 0, 0));
    IntVect ivpf0(AMREX_D_DECL(0, 0, 0));
    IntVect ivpf1(AMREX_D_DECL(0, 0, 0));
    if (has_left) {
        const int imf0 = ratio[0] * (ic - 1);
        const int imf1 = imf0 + 1;
        ivmf0 = IntVect(AMREX_D_DECL(imf0, jc, kc));
        ivmf1 = IntVect(AMREX_D_DECL(imf1, jc, kc));
    }
    if (has_right) {
        const int ipf0 = ratio[0] * (ic + 1);
        const int ipf1 = ipf0 + 1;
        ivpf0 = IntVect(AMREX_D_DECL(ipf0, jc, kc));
        ivpf1 = IntVect(AMREX_D_DECL(ipf1, jc, kc));
    }

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<GpuArray<Real,sd_order_hweno>,5> Um0_h{}, Um1_h{}, U00_h{}, U01_h{}, Up0_h{}, Up1_h{};
        for (int h = 0; h < 5; ++h) {
            const int comp = hweno_hydro_comp(h);
            for (int s = 0; s < sd_order_hweno; ++s) {
                const int off = box2point(s, line, 0, comp);
                U00_h[h][s] = srcarr(iv0f0[0], iv0f0[1], iv0f0[2], off);
                U01_h[h][s] = srcarr(iv0f1[0], iv0f1[1], iv0f1[2], off);
                if (has_left) {
                    Um0_h[h][s] = srcarr(ivmf0[0], ivmf0[1], ivmf0[2], off);
                    Um1_h[h][s] = srcarr(ivmf1[0], ivmf1[1], ivmf1[2], off);
                }
                if (has_right) {
                    Up0_h[h][s] = srcarr(ivpf0[0], ivpf0[1], ivpf0[2], off);
                    Up1_h[h][s] = srcarr(ivpf1[0], ivpf1[1], ivpf1[2], off);
                }
            }
        }

        GpuArray<RestrictInput1D,5> rin_h{};
        if (use_central_hweno) {
            for (int h = 0; h < 5; ++h) {
                const auto m_m5q4 = NodalToMoments1D(Um0_h[h].data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1_h[h].data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00_h[h].data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01_h[h].data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0_h[h].data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1_h[h].data(), hfine);
                rin_h[h] = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                                m_p1q4, m_p3q4, m_p5q4);
            }
        }

        GpuArray<GpuArray<Real,3>,5> bH_cons{};
        bool hydro_ready = false;
        if (use_central_hweno) {
            GpuArray<Real,5> Uj{};
            Uj[0] = rin_h[0].Uj;
            Uj[1] = rin_h[1].Uj;
            Uj[2] = rin_h[2].Uj;
            Uj[3] = rin_h[3].Uj;
            Uj[4] = rin_h[4].Uj;
            GpuArray<GpuArray<Real,5>,5> L{}, R{};
            GpuArray<Real,5> lambda{};
            hweno_build_eigensystem_ideal5(L, R, lambda, Uj, 0, eos_gamma);

            GpuArray<GpuArray<Real,3>,5> bH_hydro{};
            for (int a = 0; a < 5; ++a) {
                Real ujm1 = Real(0.0), uj = Real(0.0), ujp1 = Real(0.0);
                Real gjm1 = Real(0.0), gjp1 = Real(0.0);
                for (int h = 0; h < 5; ++h) {
                    ujm1 += L[a][h] * rin_h[h].Ujm1;
                    uj   += L[a][h] * rin_h[h].Uj;
                    ujp1 += L[a][h] * rin_h[h].Ujp1;
                    gjm1 += L[a][h] * rin_h[h].Gjm1_1;
                    gjp1 += L[a][h] * rin_h[h].Gjp1_1;
                }
                const auto cand = Build3Candidates(ujm1, uj, ujp1, gjm1, gjp1, hcrse);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], SmoothRegion::Parent);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH_hydro[a][m] += omega[kk] * cand[kk][m];
                    }
                }
            }
            for (int h = 0; h < 5; ++h) {
                for (int m = 0; m < 3; ++m) {
                    for (int a = 0; a < 5; ++a) {
                        bH_cons[h][m] += R[h][a] * bH_hydro[a][m];
                    }
                }
            }
            hydro_ready = true;
        }

        for (int nv = 0; nv < nvar; ++nv) {
            GpuArray<Real,3> bH{};
            if (hydro_ready && hweno_is_hydro_comp(nv)) {
                const int h = (nv == MRHO) ? 0
                            : (nv == MU)   ? 1
                            : (nv == MV)   ? 2
                            : (nv == MW)   ? 3
                                           : 4;
                bH = bH_cons[h];
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(s, line, 0, nv);
                    crsearr(i, j, k, off) = EvalQuadratic(bH, xi_sol[s]/2.0);
                }
                continue;
            }

            GpuArray<Real,sd_order_hweno> Um0{}, Um1{}, U00{}, U01{}, Up0{}, Up1{};
            for (int s = 0; s < sd_order_hweno; ++s) {
                const int off = box2point(s, line, 0, nv);
                U00[s] = srcarr(iv0f0[0], iv0f0[1], iv0f0[2], off);
                U01[s] = srcarr(iv0f1[0], iv0f1[1], iv0f1[2], off);
                if (has_left) {
                    Um0[s] = srcarr(ivmf0[0], ivmf0[1], ivmf0[2], off);
                    Um1[s] = srcarr(ivmf1[0], ivmf1[1], ivmf1[2], off);
                }
                if (has_right) {
                    Up0[s] = srcarr(ivpf0[0], ivpf0[1], ivpf0[2], off);
                    Up1[s] = srcarr(ivpf1[0], ivpf1[1], ivpf1[2], off);
                }
            }

            if (use_central_hweno) {
                const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
                const auto rin = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                                      m_p1q4, m_p3q4, m_p5q4);
                const auto cand = Build3Candidates(rin.Ujm1, rin.Uj, rin.Ujp1,
                                                   rin.Gjm1_1, rin.Gjp1_1, hcrse);
                GpuArray<Real,3> beta{};
                for (int kk = 0; kk < 3; ++kk) {
                    beta[kk] = BetaFromQuadratic(cand[kk], SmoothRegion::Parent);
                }
                const auto omega = ZWeights3(beta);
                for (int m = 0; m < 3; ++m) {
                    for (int kk = 0; kk < 3; ++kk) {
                        bH[m] += omega[kk] * cand[kk][m];
                    }
                }
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(s, line, 0, nv);
                    crsearr(i, j, k, off) = EvalQuadratic(bH, xi_sol[s]/2.0);
                }
                continue;
            }

            if (use_left_biased_p3) {
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
                const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
                const Real Uj = Real(0.5) * (m_m1q4.ubar + m_p1q4.ubar);
                const Real Ujp1 = Real(0.5) * (m_p3q4.ubar + m_p5q4.ubar);
                const Real a = Ujp1 - Uj;
                const Real b = Uj;
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(s, line, 0, nv);
                    const Real xi = xi_sol[s] / Real(2.0);
                    crsearr(i, j, k, off) = a * xi + b;
                }
            } else if (use_right_biased_p2) {
                const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
                const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
                const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
                const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
                const Real Ujm1 = Real(0.5) * (m_m5q4.ubar + m_m3q4.ubar);
                const Real Uj = Real(0.5) * (m_m1q4.ubar + m_p1q4.ubar);
                const Real a = Uj - Ujm1;
                const Real b = Uj;
                for (int s = 0; s < sd_order_hweno; ++s) {
                    const int off = box2point(s, line, 0, nv);
                    const Real xi = xi_sol[s] / Real(2.0);
                    crsearr(i, j, k, off) = a * xi + b;
                }
            } else {
                for (int s = 0; s < sd_order_hweno; ++s) {
                    Real val = Real(0.0);
                    for (int q = 0; q < sd_order_hweno; ++q) {
                        val += P1DRestrictX(0, s, q) * U00[q];
                        val += P1DRestrictX(1, s, q) * U01[q];
                    }
                    const int off = box2point(s, line, 0, nv);
                    crsearr(i, j, k, off) = val;
                }
            }
        }
    }
}

void
HermiteWENO2D::interp (const FArrayBox& crse,
                       int              crse_comp,
                       FArrayBox&       fine,
                       int              fine_comp,
                       int              ncomp,
                       const Box&       fine_region,
                       const IntVect&   ratio,
                       const Geometry&  crse_geom,
                       const Geometry&  fine_geom,
                       Vector<BCRec> const& bcr,
                       int              actual_comp,
                       int              actual_state,
                       RunOn            runon)
{
    BL_PROFILE("HermiteWENO2D::interp()");
    amrex::ignore_unused(fine_geom, bcr, actual_comp, actual_state);

    AMREX_ASSERT(ratio == 2);

    const Box target_fine_region = fine_region & fine.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);
    const int nvar = ncomp / sd_space_hweno;
    Real eos_gamma = Real(1.4);
    Real eps_rho = Real(1.0e-10);
    Real eps_p = Real(1.0e-10);

    Array4<Real const> const& carr = crse.const_array(crse_comp);
    Array4<Real>       const& farr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::coarsen(target_fine_region, IntVect(ratio[0],ratio[1],1));
    bz.grow(IntVect(1,1,0));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();

    const Real hz = crse_geom.CellSize(2);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        hweno_interp_z(i, j, k, n, tmpzarr, carr, ratio, hz);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::coarsen(target_fine_region, IntVect(AMREX_D_DECL(ratio[0],1,1)));
    by.grow(IntVect(AMREX_D_DECL(1,0,0))); // halo for x-direction stage
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = carr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif

    const Real hy = crse_geom.CellSize(1);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_interp_y(i, j, k, n, tmpyarr, srcarr, ratio, hy);
    });
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, 1, i, j, k, n,
    {
        hweno_apply_positivity_cell(tmpyarr, i, j, k, nvar, eos_gamma, eps_rho, eps_p);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = carr;
#else
    srcarr = tmpy.const_array();
#endif
    const Real hx = crse_geom.CellSize(0);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_interp_x(i, j, k, n, farr, srcarr, ratio, hx);
    });
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, 1, i, j, k, n,
    {
        hweno_apply_positivity_cell(farr, i, j, k, nvar, eos_gamma, eps_rho, eps_p);
    });
}

void
HermiteWENO2D::restrict (const FArrayBox& fine,
                         int              fine_comp,
                         FArrayBox&       crse,
                         int              crse_comp,
                         int              ncomp,
                         const Box&       crse_region,
                         const IntVect&   ratio,
                         const Geometry&  fine_geom,
                         const Geometry&  crse_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon)
{
    BL_PROFILE("HermiteWENO2D::restrict()");
    amrex::ignore_unused(bcr, actual_comp, actual_state);

    AMREX_ASSERT(ratio == 2);
    AMREX_ASSERT(ncomp % sd_space_hweno == 0);

    const Box target_crse_region = crse_region & crse.box();
    // const Box inner_box = amrex::grow(target_crse_region, -1);
    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);
    const int nvar = ncomp / sd_space_hweno;
    Real eos_gamma = Real(1.4);
    Real eps_rho = Real(1.0e-10);
    Real eps_p = Real(1.0e-10);

    Array4<Real const> const& finearr = fine.const_array(fine_comp);
    Array4<Real>       const& crsearr = crse.array(crse_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::refine(target_crse_region, IntVect(ratio[0],ratio[1],1));
    FArrayBox tmpz(bz, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpz_eli;
    if (run_on_gpu) { tmpz_eli = tmpz.elixir(); }
#endif
    Array4<Real> const& tmpzarr = tmpz.array();

    const Real hz_f = fine_geom.CellSize(2);
    const Real hz_c = crse_geom.CellSize(2);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        hweno_restrict_z(i, j, k, n, tmpzarr, finearr, ratio, hz_f, hz_c);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::refine(target_crse_region, IntVect(AMREX_D_DECL(ratio[0],1,1)));
    by.grow(IntVect(AMREX_D_DECL(ratio[0],0,0)));
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = finearr;
#else
    Array4<Real const> srcarr = tmpz.const_array(); // TODO(hweno-restrict-3d): add z-stage.
#endif

    const Real hy_f = fine_geom.CellSize(1);
    const Real hy_c = crse_geom.CellSize(1);
    const int j_valid_lo = target_crse_region.smallEnd(1);
    const int j_valid_hi = target_crse_region.bigEnd(1);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_restrict_y(i, j, k, n, tmpyarr, srcarr, ratio, hy_f, hy_c);
    });
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, 1, i, j, k, n,
    {
        hweno_apply_positivity_cell(tmpyarr, i, j, k, nvar, eos_gamma, eps_rho, eps_p);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = finearr;
#else
    srcarr = tmpy.const_array();
#endif
    const Real hx_f = fine_geom.CellSize(0);
    const Real hx_c = crse_geom.CellSize(0);
    const int i_valid_lo = target_crse_region.smallEnd(0);
    const int i_valid_hi = target_crse_region.bigEnd(0);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_restrict_x(i, j, k, n, crsearr, srcarr, ratio, hx_f, hx_c);
    });
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, 1, i, j, k, n,
    {
        hweno_apply_positivity_cell(crsearr, i, j, k, nvar, eos_gamma, eps_rho, eps_p);
    });
}

}
