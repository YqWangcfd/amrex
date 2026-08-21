
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_ParallelDescriptor.H>

#include <climits>
#include <cmath>
#include <sstream>
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
    void hweno_check_source_access(
        Array4<Real const> const& srcarr,
        IntVect const& ivm, IntVect const& ivc, IntVect const& ivp,
        int off, int direction, int i, int j, int k, int n,
        int child, int line, int solution_point) noexcept
    {
        const bool valid = srcarr.p != nullptr
            && srcarr.contains(ivm)
            && srcarr.contains(ivc)
            && srcarr.contains(ivp)
            && off >= 0 && off < srcarr.nComp();
        if (valid) {
            return;
        }

        AMREX_IF_ON_DEVICE((
            const Dim3 ivm3 = ivm.dim3();
            const Dim3 ivc3 = ivc.dim3();
            const Dim3 ivp3 = ivp.dim3();
            AMREX_DEVICE_PRINTF(
                "[HWENO_PROLONG_ACCESS_OOB] dir=%d fine=(%d,%d,%d) "
                "var=%d child=%d line=%d point=%d off=%d ncomp=%d "
                "ivm=(%d,%d,%d) ivc=(%d,%d,%d) ivp=(%d,%d,%d) "
                "src=(%d:%d,%d:%d,%d:%d)\n",
                direction, i, j, k, n, child, line, solution_point,
                off, srcarr.nComp(),
                ivm3.x, ivm3.y, ivm3.z,
                ivc3.x, ivc3.y, ivc3.z,
                ivp3.x, ivp3.y, ivp3.z,
                srcarr.begin.x, srcarr.end.x-1,
                srcarr.begin.y, srcarr.end.y-1,
                srcarr.begin.z, srcarr.end.z-1);
            amrex::Abort();
        ))
        AMREX_IF_ON_HOST((
            std::ostringstream message;
            message << "[HWENO_PROLONG_ACCESS_OOB] rank="
                    << ParallelDescriptor::MyProc()
                    << " dir=" << (direction == 0 ? "x" : "y")
                    << " fine=(" << i << "," << j << "," << k << ")"
                    << " var=" << n << " child=" << child
                    << " line=" << line << " point=" << solution_point
                    << " off=" << off << " ncomp=" << srcarr.nComp()
                    << " ivm=" << ivm << " ivc=" << ivc << " ivp=" << ivp
                    << " src=(" << srcarr.begin.x << ":" << srcarr.end.x-1
                    << "," << srcarr.begin.y << ":" << srcarr.end.y-1
                    << "," << srcarr.begin.z << ":" << srcarr.end.z-1
                    << ")";
            amrex::Abort(message.str());
        ))
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_pp_internal_energy (Real rho, Real mx, Real my, Real mz, Real e) noexcept
    {
        if (!(rho > Real(0.0)) || !amrex::Math::isfinite(rho)) {
            return -AMREX_REAL_MAX;
        }
        return e - Real(0.5) * (mx*mx + my*my + mz*mz) / rho;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_pp_lagrange_basis (int node, Real x) noexcept
    {
        Real value = Real(1.0);
        for (int other = 0; other < sd_ORDER; ++other) {
            if (other != node) {
                value *= (x - xi_sol[other]) / (xi_sol[node] - xi_sol[other]);
            }
        }
        return value;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_pp_sample_child (Array4<Real const> const& fine, int i, int j, int k,
                              int var, int point) noexcept
    {
        const int base = var * sd_SPACE;
        if (point < sd_SPACE) {
            return fine(i,j,k,base + point);
        }

        point -= sd_SPACE;
        if (point < sd_edge_ORDER * sd_ORDER) {
            const int iy = point / sd_edge_ORDER;
            const int fx = point - iy * sd_edge_ORDER;
            Real value = Real(0.0);
            for (int ix = 0; ix < sd_ORDER; ++ix) {
                value += amr_pp_lagrange_basis(ix, xi_flx[fx])
                    * fine(i,j,k,base + box2point(ix,iy,0,0));
            }
            return value;
        }

        point -= sd_edge_ORDER * sd_ORDER;
        const int ix = point / sd_edge_ORDER;
        const int fy = point - ix * sd_edge_ORDER;
        Real value = Real(0.0);
        for (int iy = 0; iy < sd_ORDER; ++iy) {
            value += amr_pp_lagrange_basis(iy, xi_flx[fy])
                * fine(i,j,k,base + box2point(ix,iy,0,0));
        }
        return value;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_pp_parent_average (Array4<Real const> const& crse,
                                int i, int j, int k, int var) noexcept
    {
        Real value = Real(0.0);
        const int base = var * sd_SPACE;
        for (int iy = 0; iy < sd_ORDER; ++iy) {
            for (int ix = 0; ix < sd_ORDER; ++ix) {
                value += SDwgh1D[ix] * SDwgh1D[iy]
                    * crse(i,j,k,base + box2point(ix,iy,0,0));
            }
        }
        return value;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_pp_child_average (Array4<Real const> const& child,
                               int i, int j, int k, int var,
                               IntVect const& child_ratio) noexcept
    {
        Real value = Real(0.0);
        const Real weight = Real(1.0)
            / Real(child_ratio[0] * child_ratio[1]);
        for (int joff = 0; joff < child_ratio[1]; ++joff) {
            for (int ioff = 0; ioff < child_ratio[0]; ++ioff) {
                value += weight * amr_pp_parent_average(
                    child,
                    i * child_ratio[0] + ioff,
                    j * child_ratio[1] + joff,
                    k, var);
            }
        }
        return value;
    }

    using AMRRestrictionPPState = GpuArray<Real, NEQ>;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real amr_restriction_pressure (AMRRestrictionPPState const& state,
                                   Real gamma) noexcept
    {
        return (gamma - Real(1.0)) * amr_pp_internal_energy(
            state[MRHO], state[MU], state[MV], state[MW], state[ME]);
    }

    AMRRestrictionPPCounters apply_restriction_positivity (
        FArrayBox& crse, Box const& region, int crse_comp, int nvar,
        FArrayBox const* source, int source_comp,
        IntVect const& source_ratio,
        bool enabled, Real gamma, Real eps_rho, Real eps_p, RunOn runon)
    {
        AMRRestrictionPPCounters counts;
        if (!enabled || !region.ok()) {
            return counts;
        }

#if (AMREX_SPACEDIM != 2)
        amrex::ignore_unused(
            crse, crse_comp, nvar, source, source_comp, source_ratio,
            gamma, eps_rho, eps_p, runon);
        amrex::Abort(
            "AMR restriction PP for SD-packed states is implemented only in 2D.");
        return counts;
#else
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nvar == NEQ,
            "AMR restriction PP requires a complete SD-packed state.");

        IArrayBox flags(region, 4);
        flags.setVal(0);
        auto const& state = crse.array(crse_comp);
        auto const& flag = flags.array();
        Array4<Real const> source_state;
        const bool has_source = source != nullptr;
        if (has_source) {
            source_state = source->const_array(source_comp);
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, region, i, j, k,
        {
            AMRRestrictionPPState Ubar{};
            bool average_nonfinite = false;
            bool raw_nonfinite = false;
            Real rho_min = AMREX_REAL_MAX;

            if (has_source) {
                for (int n = 0; n < NEQ; ++n) {
                    Ubar[n] = amr_pp_child_average(
                        source_state, i, j, k, n, source_ratio);
                    average_nonfinite = average_nonfinite
                        || !amrex::Math::isfinite(Ubar[n]);
                }
            }

            for (int iy = 0; iy < sd_ORDER; ++iy) {
                for (int ix = 0; ix < sd_ORDER; ++ix) {
                    const int off = box2point(ix,iy,0,0);
                    const Real weight = SDwgh1D[ix] * SDwgh1D[iy];
                    for (int n = 0; n < NEQ; ++n) {
                        const Real value = state(i,j,k,n*sd_SPACE + off);
                        raw_nonfinite = raw_nonfinite
                            || !amrex::Math::isfinite(value);
                        if (!has_source) {
                            Ubar[n] += weight * value;
                        }
                    }
                    rho_min = amrex::min(
                        rho_min, state(i,j,k,MRHO*sd_SPACE + off));
                }
            }

            for (int n = 0; n < NEQ; ++n) {
                average_nonfinite = average_nonfinite
                    || !amrex::Math::isfinite(Ubar[n]);
            }
            bool species_average_invalid = false;
            for (int n = 0; n < NSP; ++n) {
                species_average_invalid = species_average_invalid
                    || Ubar[n] < Real(0.0);
            }
            const Real rho_bar = Ubar[MRHO];
            const Real p_bar = amr_restriction_pressure(Ubar, gamma);
            if (average_nonfinite || species_average_invalid
                || !amrex::Math::isfinite(p_bar)
                || rho_bar < eps_rho || p_bar < eps_p) {
                flag(i,j,k,2) = 1;
                return;
            }

            Real theta_rho = Real(1.0);
            if (raw_nonfinite) {
                theta_rho = Real(0.0);
            } else if (rho_min < eps_rho) {
                const Real denominator = rho_bar - rho_min;
                theta_rho = denominator > Real(1.0e-30)
                    ? amrex::min(
                        Real(1.0), (rho_bar - eps_rho) / denominator)
                    : Real(0.0);
            }

            Real theta_species = Real(1.0);
            if (!raw_nonfinite) {
                for (int off = 0; off < sd_SPACE; ++off) {
                    for (int n = 0; n < NSP; ++n) {
                        const Real raw = state(i,j,k,n*sd_SPACE + off);
                        const Real density_limited = Ubar[n]
                            + theta_rho * (raw - Ubar[n]);
                        if (density_limited < Real(0.0)) {
                            const Real denominator = Ubar[n] - density_limited;
                            theta_species = denominator > Real(1.0e-30)
                                ? amrex::min(
                                    theta_species, Ubar[n] / denominator)
                                : Real(0.0);
                        }
                    }
                }
            }

            Real theta_pressure = Real(1.0);
            if (raw_nonfinite) {
                theta_pressure = Real(0.0);
            }
            for (int off = 0; off < sd_SPACE && !raw_nonfinite; ++off) {
                AMRRestrictionPPState Uhat{};
                for (int n = 0; n < NEQ; ++n) {
                    const Real raw = state(i,j,k,n*sd_SPACE + off);
                    const Real density_limited = Ubar[n]
                        + theta_rho * (raw - Ubar[n]);
                    Uhat[n] = Ubar[n]
                        + theta_species * (density_limited - Ubar[n]);
                }

                const Real p_hat = amr_restriction_pressure(Uhat, gamma);
                if (!amrex::Math::isfinite(p_hat)
                    || Uhat[MRHO] < eps_rho || p_hat < eps_p) {
                    Real lo = Real(0.0);
                    Real hi = Real(1.0);
                    for (int iter = 0; iter < 64; ++iter) {
                        const Real mid = Real(0.5) * (lo + hi);
                        AMRRestrictionPPState Umid{};
                        for (int n = 0; n < NEQ; ++n) {
                            Umid[n] = Ubar[n] + mid * (Uhat[n] - Ubar[n]);
                        }
                        const Real p_mid = amr_restriction_pressure(Umid, gamma);
                        if (amrex::Math::isfinite(p_mid)
                            && Umid[MRHO] >= eps_rho && p_mid >= eps_p) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    theta_pressure = amrex::min(theta_pressure, lo);
                }
            }

            flag(i,j,k,0) = !raw_nonfinite && theta_rho < Real(1.0);
            flag(i,j,k,1) = !raw_nonfinite && theta_pressure < Real(1.0);
            flag(i,j,k,3) = flag(i,j,k,0) || flag(i,j,k,1)
                || theta_species < Real(1.0);
            flag(i,j,k,3) = flag(i,j,k,3) || raw_nonfinite;

            if (flag(i,j,k,3) != 0) {
                for (int off = 0; off < sd_SPACE; ++off) {
                    for (int n = 0; n < NEQ; ++n) {
                        if (raw_nonfinite) {
                            state(i,j,k,n*sd_SPACE + off) = Ubar[n];
                            continue;
                        }
                        const Real raw = state(i,j,k,n*sd_SPACE + off);
                        const Real density_limited = Ubar[n]
                            + theta_rho * (raw - Ubar[n]);
                        const Real species_limited = Ubar[n]
                            + theta_species * (density_limited - Ubar[n]);
                        state(i,j,k,n*sd_SPACE + off) = Ubar[n]
                            + theta_pressure * (species_limited - Ubar[n]);
                    }
                }
            }
        });

        const int invalid = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(region, 2)
            : flags.sum<RunOn::Host>(region, 2);
        if (invalid != 0) {
            amrex::Abort(
                "AMR restriction PP requires finite point values and an admissible coarse-cell average.");
        }

        counts.rho = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(region, 0)
            : flags.sum<RunOn::Host>(region, 0);
        counts.pressure = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(region, 1)
            : flags.sum<RunOn::Host>(region, 1);
        counts.any = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(region, 3)
            : flags.sum<RunOn::Host>(region, 3);
        return counts;
#endif
    }

    constexpr int amr_pp_check_points = sd_SPACE + 2 * sd_edge_ORDER * sd_ORDER;

    AMRProlongationPPCounters apply_parentwise_prolongation_pp (
        FArrayBox& raw_fine, FArrayBox const& crse, Box const& parent_region,
        int crse_comp, IntVect const& child_ratio, int nvar,
        bool limit_positivity, bool check_flux_points,
        Real gamma, Real eps_rho, Real eps_p,
        IArrayBox const* forced_nonfinite, RunOn runon)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !limit_positivity || nvar == NEQ,
            "AMR prolongation PP requires a complete SD-packed state.");
        IArrayBox flags(parent_region, 5);
        flags.setVal(0);
        auto const& fine = raw_fine.array();
        auto const& fine_const = raw_fine.const_array();
        auto const& coarse = crse.const_array(crse_comp);
        auto const& flag = flags.array();
        Array4<int const> forced;
        const bool has_forced = forced_nonfinite != nullptr;
        if (has_forced) {
            forced = forced_nonfinite->const_array();
        }
        const Real eps_internal = eps_p / (gamma - Real(1.0));
        const int number_of_check_points = check_flux_points
            ? amr_pp_check_points : sd_SPACE;

        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, parent_region, i, j, k,
        {
            AMRRestrictionPPState Ubar{};
            bool parent_average_nonfinite = false;
            for (int n = 0; n < nvar; ++n) {
                const Real average = amr_pp_parent_average(
                    coarse, i, j, k, n);
                parent_average_nonfinite = parent_average_nonfinite
                    || !amrex::Math::isfinite(average);
                if (limit_positivity) {
                    Ubar[n] = average;
                }
            }

            Real rho_bar = Real(0.0);
            Real internal_bar = Real(0.0);
            if (limit_positivity) {
                rho_bar = Ubar[MRHO];
                const Real mx_bar = Ubar[MU];
                const Real my_bar = Ubar[MV];
                const Real mz_bar = Ubar[MW];
                const Real e_bar = Ubar[ME];
                internal_bar = amr_pp_internal_energy(
                    rho_bar, mx_bar, my_bar, mz_bar, e_bar);
            }

            if (parent_average_nonfinite
                || (limit_positivity
                    && (!amrex::Math::isfinite(rho_bar)
                        || !amrex::Math::isfinite(internal_bar)
                        || rho_bar < eps_rho || internal_bar < eps_internal))) {
                flag(i,j,k,3) = 1;
                return;
            }
            if (limit_positivity) {
                for (int n = 0; n < NSP; ++n) {
                    if (Ubar[n] < Real(0.0)) {
                        flag(i,j,k,3) = 1;
                        return;
                    }
                }
            }

            const bool forced_fallback = has_forced && forced(i,j,k,0) != 0;
            bool raw_nonfinite = false;
            Real rho_min = AMREX_REAL_MAX;
            for (int joff = 0; joff < child_ratio[1]; ++joff) {
                for (int ioff = 0; ioff < child_ratio[0]; ++ioff) {
                    const int fi = i * child_ratio[0] + ioff;
                    const int fj = j * child_ratio[1] + joff;
                    for (int point = 0; point < number_of_check_points; ++point) {
                        for (int n = 0; n < nvar; ++n) {
                            raw_nonfinite = raw_nonfinite
                                || !amrex::Math::isfinite(
                                amr_pp_sample_child(fine_const, fi, fj, k, n, point));
                        }
                        if (limit_positivity) {
                            rho_min = amrex::min(
                                rho_min,
                                amr_pp_sample_child(
                                    fine_const, fi, fj, k, MRHO, point));
                        }
                    }
                }
            }
            bool nonfinite = forced_fallback || raw_nonfinite;

            Real theta_rho = Real(1.0);
            if (nonfinite) {
                theta_rho = Real(0.0);
            } else if (limit_positivity && rho_min < eps_rho) {
                const Real denominator = rho_bar - rho_min;
                theta_rho = denominator > Real(1.0e-30)
                    ? amrex::min(Real(1.0), (rho_bar - eps_rho) / denominator)
                    : Real(0.0);
            }

            Real theta_species = Real(1.0);
            if (limit_positivity && !nonfinite) {
                for (int joff = 0; joff < child_ratio[1]; ++joff) {
                    for (int ioff = 0; ioff < child_ratio[0]; ++ioff) {
                        const int fi = i * child_ratio[0] + ioff;
                        const int fj = j * child_ratio[1] + joff;
                        for (int point = 0; point < number_of_check_points; ++point) {
                            for (int n = 0; n < NSP; ++n) {
                                const Real raw = amr_pp_sample_child(
                                    fine_const, fi, fj, k, n, point);
                                const Real density_limited = Ubar[n]
                                    + theta_rho * (raw - Ubar[n]);
                                if (density_limited < Real(0.0)) {
                                    const Real denominator =
                                        Ubar[n] - density_limited;
                                    theta_species = denominator > Real(1.0e-30)
                                        ? amrex::min(
                                            theta_species,
                                            Ubar[n] / denominator)
                                        : Real(0.0);
                                }
                            }
                        }
                    }
                }
            }

            const Real theta_admissible = theta_rho * theta_species;

            Real internal_min = AMREX_REAL_MAX;
            if (limit_positivity && !nonfinite) {
                for (int joff = 0; joff < child_ratio[1]; ++joff) {
                    for (int ioff = 0; ioff < child_ratio[0]; ++ioff) {
                        const int fi = i * child_ratio[0] + ioff;
                        const int fj = j * child_ratio[1] + joff;
                        for (int point = 0; point < number_of_check_points; ++point) {
                            const Real rho_raw = amr_pp_sample_child(
                                fine_const, fi, fj, k, MRHO, point);
                            const Real rho = rho_bar
                                + theta_admissible * (rho_raw - rho_bar);
                            const Real mx = Ubar[MU] + theta_admissible * (
                                amr_pp_sample_child(
                                    fine_const, fi, fj, k, MU, point) - Ubar[MU]);
                            const Real my = Ubar[MV] + theta_admissible * (
                                amr_pp_sample_child(
                                    fine_const, fi, fj, k, MV, point) - Ubar[MV]);
                            const Real mz = Ubar[MW] + theta_admissible * (
                                amr_pp_sample_child(
                                    fine_const, fi, fj, k, MW, point) - Ubar[MW]);
                            const Real energy = Ubar[ME] + theta_admissible * (
                                amr_pp_sample_child(
                                    fine_const, fi, fj, k, ME, point) - Ubar[ME]);
                            const Real internal = amr_pp_internal_energy(
                                rho, mx, my, mz, energy);
                            if (!amrex::Math::isfinite(internal)) {
                                raw_nonfinite = true;
                                nonfinite = true;
                            }
                            internal_min = amrex::min(internal_min, internal);
                        }
                    }
                }
            }

            Real theta_internal = Real(1.0);
            if (nonfinite) {
                theta_rho = Real(0.0);
                theta_internal = Real(0.0);
            } else if (limit_positivity && internal_min < eps_internal) {
                const Real denominator = internal_bar - internal_min;
                theta_internal = denominator > Real(1.0e-30)
                    ? amrex::min(Real(1.0),
                        (internal_bar - eps_internal) / denominator)
                    : Real(0.0);
            }

            flag(i,j,k,0) = limit_positivity && !nonfinite
                && theta_rho < Real(1.0);
            flag(i,j,k,1) = limit_positivity && !nonfinite
                && theta_internal < Real(1.0);
            flag(i,j,k,2) = raw_nonfinite;
            flag(i,j,k,4) = flag(i,j,k,0) || flag(i,j,k,1)
                || theta_species < Real(1.0);

            if (theta_admissible < Real(1.0)
                || theta_internal < Real(1.0)) {
                for (int joff = 0; joff < child_ratio[1]; ++joff) {
                    for (int ioff = 0; ioff < child_ratio[0]; ++ioff) {
                        const int fi = i * child_ratio[0] + ioff;
                        const int fj = j * child_ratio[1] + joff;
                        for (int n = 0; n < nvar; ++n) {
                            const Real average = limit_positivity
                                ? Ubar[n]
                                : amr_pp_parent_average(coarse, i, j, k, n);
                            for (int off = 0; off < sd_SPACE; ++off) {
                                const int comp = n * sd_SPACE + off;
                                if (nonfinite) {
                                    fine(fi,fj,k,comp) = average;
                                    continue;
                                }
                                const Real value = fine(fi,fj,k,comp);
                                const Real density_limited = average
                                    + theta_admissible * (value - average);
                                fine(fi,fj,k,comp) = average
                                    + theta_internal * (density_limited - average);
                            }
                        }
                    }
                }
            }
        });

        int invalid_parent = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(parent_region, 3)
            : flags.sum<RunOn::Host>(parent_region, 3);
        if (invalid_parent != 0) {
            amrex::Abort(
                "AMR prolongation requires finite parent averages and, when PP is enabled, an admissible gas-state average.");
        }

        AMRProlongationPPCounters counts;
        counts.rho = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(parent_region, 0)
            : flags.sum<RunOn::Host>(parent_region, 0);
        counts.pressure = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(parent_region, 1)
            : flags.sum<RunOn::Host>(parent_region, 1);
        counts.any = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(parent_region, 4)
            : flags.sum<RunOn::Host>(parent_region, 4);
        counts.nonfinite = runon == RunOn::Gpu
            ? flags.sum<RunOn::Device>(parent_region, 2)
            : flags.sum<RunOn::Host>(parent_region, 2);
        return counts;
    }

    Long hweno_check_intermediate_finite (
        FArrayBox& intermediate, FArrayBox const& crse,
        Box const& parent_region, IntVect const& ratio,
        int crse_comp, int nvar, IArrayBox& nonfinite_parent, RunOn runon)
    {
        auto const& tmp = intermediate.array();
        auto const& coarse = crse.const_array(crse_comp);
        auto const& bad = nonfinite_parent.array();

        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, parent_region, i, j, k,
        {
            bool nonfinite = false;
            for (int joff = 0; joff < ratio[1]; ++joff) {
                const int fj = j * ratio[1] + joff;
                for (int n = 0; n < nvar; ++n) {
                    for (int off = 0; off < sd_SPACE; ++off) {
                        nonfinite = nonfinite
                            || !amrex::Math::isfinite(
                                tmp(i,fj,k,n*sd_SPACE + off));
                    }
                }
            }
            bad(i,j,k,0) = nonfinite;
            if (nonfinite) {
                for (int joff = 0; joff < ratio[1]; ++joff) {
                    const int fj = j * ratio[1] + joff;
                    for (int n = 0; n < nvar; ++n) {
                        const Real average = amr_pp_parent_average(
                            coarse, i, j, k, n);
                        for (int off = 0; off < sd_SPACE; ++off) {
                            tmp(i,fj,k,n*sd_SPACE + off) = average;
                        }
                    }
                }
            }
        });
        return runon == RunOn::Gpu
            ? nonfinite_parent.sum<RunOn::Device>(parent_region, 0)
            : nonfinite_parent.sum<RunOn::Host>(parent_region, 0);
    }

    void amr_assert_fab_finite (FArrayBox const& fab, Box const& region,
                                int comp, int ncomp, RunOn runon,
                                char const* message)
    {
        IArrayBox bad(region, 1);
        bad.setVal(0);
        auto const& src = fab.const_array(comp);
        auto const& flag = bad.array();
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, region, i, j, k,
        {
            bool nonfinite = false;
            for (int n = 0; n < ncomp; ++n) {
                nonfinite = nonfinite || !amrex::Math::isfinite(src(i,j,k,n));
            }
            flag(i,j,k,0) = nonfinite;
        });
        const int count = runon == RunOn::Gpu
            ? bad.sum<RunOn::Device>(region, 0)
            : bad.sum<RunOn::Host>(region, 0);
        if (count != 0) {
            amrex::Abort(message);
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
Mortar2D::configure_amr_transfer_pp (bool enabled, Real gamma,
                                    Real eps_rho, Real eps_p) noexcept
{
    m_prolongation_pp_enabled = enabled;
    m_pp_gamma = gamma;
    m_pp_eps_rho = eps_rho;
    m_pp_eps_p = eps_p;
}

AMRProlongationPPCounters
Mortar2D::take_prolongation_pp_counters () noexcept
{
    return {
        m_pp_rho_events.exchange(0),
        m_pp_pressure_events.exchange(0),
        m_pp_any_events.exchange(0),
        m_pp_nonfinite_events.exchange(0)
    };
}

AMRRestrictionPPCounters
Mortar2D::take_restriction_pp_counters () noexcept
{
    return {
        m_restrict_rho_events.exchange(0),
        m_restrict_pressure_events.exchange(0),
        m_restrict_any_events.exchange(0)
    };
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
    if (!target_fine_region.ok()) {
        return;
    }

    // amrex::Print() << "target_fine_region=" << target_fine_region   << std::endl;

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    // amrex::Print() << "crse_comp= " << crse_comp
    //                << "\nfine_comp= " << fine_comp << std::endl;
    // amrex::Abort("aa");

    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ncomp % sd_space == 0,
        "Mortar2D interpolation requires complete SD-packed variables.");
    if (type == Type::ScaleRef && m_prolongation_pp_enabled) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ncomp / sd_space == NEQ,
            "Mortar2D prolongation PP requires a complete SD-packed state.");
    }

    if (type != Type::ScaleRef || !m_prolongation_pp_enabled) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
            runon, target_fine_region, ncomp/sd_space, i, j, k, n,
        {
            mortar_interp(i,j,k,n,finearr,crsearr,ratio);
        });
        return;
    }

    const Box parent_region = amrex::coarsen(target_fine_region, ratio);
    const Box full_fine_region = amrex::refine(parent_region, ratio);
    FArrayBox raw_fine(full_fine_region, ncomp);
#ifdef AMREX_USE_GPU
    Elixir raw_fine_eli;
    if (run_on_gpu) { raw_fine_eli = raw_fine.elixir(); }
#endif
    auto const& raw = raw_fine.array();

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, full_fine_region, ncomp/sd_space, i, j, k, n,
    {
        mortar_interp(i,j,k,n,raw,crsearr,ratio);
    });

    const auto counts = apply_parentwise_prolongation_pp(
        raw_fine, crse, parent_region, crse_comp, ratio, ncomp/sd_space,
        true, true,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, nullptr, runon);
    m_pp_rho_events.fetch_add(counts.rho);
    m_pp_pressure_events.fetch_add(counts.pressure);
    m_pp_any_events.fetch_add(counts.any);
    m_pp_nonfinite_events.fetch_add(counts.nonfinite);

    auto const& raw_const = raw_fine.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, target_fine_region, ncomp, i, j, k, n,
    {
        finearr(i,j,k,n) = raw_const(i,j,k,n);
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
    if (m_prolongation_pp_enabled) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ncomp / sd_space == NEQ,
            "Mortar2D restriction PP requires a complete SD-packed state.");
    }
    auto const& destarr = crsearr;
    auto const& srcarr  = finearr;
    const Type restriction_type = type;
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_crse_region, ncomp/sd_space, i, j, k, n,
    {
        Mortar2D::mortar_restrict_with_type(
            i, j, k, n, destarr, srcarr, ratio, restriction_type);
        // mortar_restrict_flatten(i,j,k,n,destarr,srcarr,ratio);
    });
    const auto counts = apply_restriction_positivity(
        crse, target_crse_region, crse_comp, ncomp/sd_space,
        nullptr, 0, IntVect(1),
        m_prolongation_pp_enabled,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, runon);
    m_restrict_rho_events.fetch_add(counts.rho);
    m_restrict_pressure_events.fetch_add(counts.pressure);
    m_restrict_any_events.fetch_add(counts.any);
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Mortar2D::mortar_restrict(const int i, const int j, const int k, const int n,
                Array4<Real> const crsearr, Array4<const Real> const& finearr,
                const IntVect& ratio)
{
    mortar_restrict_with_type(
        i, j, k, n, crsearr, finearr, ratio, type);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Mortar2D::mortar_restrict_with_type(
                const int i, const int j, const int k, const int n,
                Array4<Real> const crsearr, Array4<const Real> const& finearr,
                const IntVect& ratio, Type restriction_type)
{
    AMREX_ASSERT(ratio[0]==2 && ratio[1]==2);
    if (restriction_type == Type::OrderRef)
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

    } else if (restriction_type == Type::ScaleRef)
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
                               ProlongWeightMode mode,
                               Real hdir) noexcept
{
    // In d-by-d interp(), y-stage loops over `by`, where x-index is already coarse.
    const int ic = i;
    const int jc = amrex::coarsen(j, ratio[1]);
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int child = j - ratio[1]*jc;
    const IntVect ivm(AMREX_D_DECL(ic, jc-1, kc));
    const IntVect ivc(AMREX_D_DECL(ic, jc,   kc));
    const IntVect ivp(AMREX_D_DECL(ic, jc+1, kc));

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(line, s, 0, n); // y-dir: line is x-node
            hweno_check_source_access(
                srcarr, ivm, ivc, ivp, off, 1,
                i, j, k, n, child, line, s);
            Um[s] = srcarr(ivm, off);
            U0[s] = srcarr(ivc, off);
            Up[s] = srcarr(ivp, off);
        }

        const auto mm = NodalToMoments1D(Um.data(), hdir);
        const auto m0 = NodalToMoments1D(U0.data(), hdir);
        const auto mp = NodalToMoments1D(Up.data(), hdir);
        const auto cand = Build4Candidates(mm.ubar, m0.ubar, mp.ubar,
                                           mm.g1, mp.g1, mm.g2, mp.g2, hdir);

        const auto pair = BuildProlongPolynomialPair(cand, m0.ubar, mode);
        const auto& polynomial = (child == 0) ? pair.left : pair.right;
        const Real ubar_parent = Real(0.5) * (
            ProlongChildAverage(pair.left, 0)
            + ProlongChildAverage(pair.right, 1));
        const Real cons_abs = std::abs(ubar_parent - m0.ubar);
        const Real cons_rel = cons_abs / amrex::max(std::abs(m0.ubar), Real(1.0e-14));
#if !defined(AMREX_USE_GPU)
        if (cons_rel > Real(1.0e-10) && ParallelDescriptor::IOProcessor()) {
            static int warn_count_y = 0;
            if (warn_count_y < 12) {
                amrex::Print() << "[HermiteWENO2D::hweno_interp_y] parent poly avg mismatch: "
                               << "abs=" << cons_abs << ", rel=" << cons_rel
                               << ", coarse=(" << ic << "," << jc << "), child=" << child
                               << ", line=" << line << ", var=" << n << "\n";
                ++warn_count_y;
            }
        }
#endif

        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(line, s, 0, n);
            tmparr(i, j, k, off) = EvalCubic(
                polynomial, hweno_child_coordinate(child, s));
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
                               ProlongWeightMode mode,
                               Real hdir) noexcept
{
    const int ic = amrex::coarsen(i, ratio[0]);
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int child = i - ratio[0]*ic;
    const IntVect ivm(AMREX_D_DECL(ic-1, jc, kc));
    const IntVect ivc(AMREX_D_DECL(ic,   jc, kc));
    const IntVect ivp(AMREX_D_DECL(ic+1, jc, kc));

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um{}, U0{}, Up{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(s, line, 0, n); // x-dir: line is y-node
            hweno_check_source_access(
                srcarr, ivm, ivc, ivp, off, 0,
                i, j, k, n, child, line, s);
            Um[s] = srcarr(ivm, off);
            U0[s] = srcarr(ivc, off);
            Up[s] = srcarr(ivp, off);
        }

        const auto mm = NodalToMoments1D(Um.data(), hdir);
        const auto m0 = NodalToMoments1D(U0.data(), hdir);
        const auto mp = NodalToMoments1D(Up.data(), hdir);
        const auto cand = Build4Candidates(mm.ubar, m0.ubar, mp.ubar,
                                           mm.g1, mp.g1, mm.g2, mp.g2, hdir);

        const auto pair = BuildProlongPolynomialPair(cand, m0.ubar, mode);
        const auto& polynomial = (child == 0) ? pair.left : pair.right;
        const Real ubar_parent = Real(0.5) * (
            ProlongChildAverage(pair.left, 0)
            + ProlongChildAverage(pair.right, 1));
        const Real cons_abs = std::abs(ubar_parent - m0.ubar);
        const Real cons_rel = cons_abs / amrex::max(std::abs(m0.ubar), Real(1.0e-14));
#if !defined(AMREX_USE_GPU)
        if (cons_rel > Real(1.0e-10) && ParallelDescriptor::IOProcessor()) {
            static int warn_count_x = 0;
            if (warn_count_x < 12) {
                amrex::Print() << "[HermiteWENO2D::hweno_interp_x] parent poly avg mismatch: "
                               << "abs=" << cons_abs << ", rel=" << cons_rel
                               << ", coarse=(" << ic << "," << jc << "), child=" << child
                               << ", line=" << line << ", var=" << n << "\n";
                ++warn_count_x;
            }
        }
#endif

        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(s, line, 0, n);
            finearr(i, j, k, off) = EvalCubic(
                polynomial, hweno_child_coordinate(child, s));
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
                                 Real hcrse) noexcept
{
    const int ic = i;
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int jmf0 = ratio[1] * (jc - 1);
    const int jmf1 = jmf0 + 1;
    const int j0f0 = ratio[1] * jc;
    const int j0f1 = j0f0 + 1;
    const int jpf0 = ratio[1] * (jc + 1);
    const int jpf1 = jpf0 + 1;

    const IntVect ivmf0(AMREX_D_DECL(ic, jmf0, kc));
    const IntVect ivmf1(AMREX_D_DECL(ic, jmf1, kc));
    const IntVect iv0f0(AMREX_D_DECL(ic, j0f0, kc));
    const IntVect iv0f1(AMREX_D_DECL(ic, j0f1, kc));
    const IntVect ivpf0(AMREX_D_DECL(ic, jpf0, kc));
    const IntVect ivpf1(AMREX_D_DECL(ic, jpf1, kc));

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um0{}, Um1{}, U00{}, U01{}, Up0{}, Up1{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(line, s, 0, n);
            Um0[s] = srcarr(ivmf0, off);
            Um1[s] = srcarr(ivmf1, off);
            U00[s] = srcarr(iv0f0, off);
            U01[s] = srcarr(iv0f1, off);
            Up0[s] = srcarr(ivpf0, off);
            Up1[s] = srcarr(ivpf1, off);
        }

        const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
        const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
        const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
        const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
        const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
        const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
        const auto rin = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                              m_p1q4, m_p3q4, m_p5q4);

        const auto cand = Build4Candidates(rin.Ujm1, rin.Uj, rin.Ujp1,
                                           rin.Gjm1_1, rin.Gjp1_1,
                                           rin.Gjm1_2, rin.Gjp1_2, hcrse);

        GpuArray<Real,4> beta{};
        for (int kk = 0; kk < 4; ++kk) {
            beta[kk] = BetaFromCubic(cand[kk], SmoothRegion::Parent);
        }
        const auto omega = JSWeights(beta);

        GpuArray<Real,4> bH{};
        for (int m = 0; m < 4; ++m) {
            for (int kk = 0; kk < 4; ++kk) {
                bH[m] += omega[kk] * cand[kk][m];
            }
        }

        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(line, s, 0, n);
            tmparr(i, j, k, off) = EvalCubic(bH, xi_sol[s]/2.0);
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
                                 Real hcrse) noexcept
{
    const int ic = i;
    const int jc = j;
    const int kc = (AMREX_SPACEDIM > 2) ? k : 0;

    const int imf0 = ratio[0] * (ic - 1);
    const int imf1 = imf0 + 1;
    const int i0f0 = ratio[0] * ic;
    const int i0f1 = i0f0 + 1;
    const int ipf0 = ratio[0] * (ic + 1);
    const int ipf1 = ipf0 + 1;

    const IntVect ivmf0(AMREX_D_DECL(imf0, jc, kc));
    const IntVect ivmf1(AMREX_D_DECL(imf1, jc, kc));
    const IntVect iv0f0(AMREX_D_DECL(i0f0, jc, kc));
    const IntVect iv0f1(AMREX_D_DECL(i0f1, jc, kc));
    const IntVect ivpf0(AMREX_D_DECL(ipf0, jc, kc));
    const IntVect ivpf1(AMREX_D_DECL(ipf1, jc, kc));

    for (int line = 0; line < sd_order_hweno; ++line) {
        GpuArray<Real,sd_order_hweno> Um0{}, Um1{}, U00{}, U01{}, Up0{}, Up1{};
        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(s, line, 0, n);
            Um0[s] = srcarr(ivmf0, off);
            Um1[s] = srcarr(ivmf1, off);
            U00[s] = srcarr(iv0f0, off);
            U01[s] = srcarr(iv0f1, off);
            Up0[s] = srcarr(ivpf0, off);
            Up1[s] = srcarr(ivpf1, off);
        }

        const auto m_m5q4 = NodalToMoments1D(Um0.data(), hfine);
        const auto m_m3q4 = NodalToMoments1D(Um1.data(), hfine);
        const auto m_m1q4 = NodalToMoments1D(U00.data(), hfine);
        const auto m_p1q4 = NodalToMoments1D(U01.data(), hfine);
        const auto m_p3q4 = NodalToMoments1D(Up0.data(), hfine);
        const auto m_p5q4 = NodalToMoments1D(Up1.data(), hfine);
        const auto rin = BuildRestrictInput1D(m_m5q4, m_m3q4, m_m1q4,
                                              m_p1q4, m_p3q4, m_p5q4);

        const auto cand = Build4Candidates(rin.Ujm1, rin.Uj, rin.Ujp1,
                                           rin.Gjm1_1, rin.Gjp1_1,
                                           rin.Gjm1_2, rin.Gjp1_2, hcrse);

        GpuArray<Real,4> beta{};
        for (int kk = 0; kk < 4; ++kk) {
            beta[kk] = BetaFromCubic(cand[kk], SmoothRegion::Parent);
        }
        const auto omega = JSWeights(beta);

        GpuArray<Real,4> bH{};
        for (int m = 0; m < 4; ++m) {
            for (int kk = 0; kk < 4; ++kk) {
                bH[m] += omega[kk] * cand[kk][m];
            }
        }

        for (int s = 0; s < sd_order_hweno; ++s) {
            const int off = box2point(s, line, 0, n);
            crsearr(i, j, k, off) = EvalCubic(bH, xi_sol[s]/2.0);
        }
    }
}

void
HermiteWENO2D::configure_amr_transfer_pp (bool enabled, Real gamma,
                                         Real eps_rho, Real eps_p) noexcept
{
    m_prolongation_pp_enabled = enabled;
    m_pp_gamma = gamma;
    m_pp_eps_rho = eps_rho;
    m_pp_eps_p = eps_p;
}

void
HermiteWENO2D::configure_prolong_weight_mode (ProlongWeightMode mode) noexcept
{
    m_prolong_weight_mode = mode;
}

AMRProlongationPPCounters
HermiteWENO2D::take_prolongation_pp_counters () noexcept
{
    return {
        m_pp_rho_events.exchange(0),
        m_pp_pressure_events.exchange(0),
        m_pp_any_events.exchange(0),
        m_pp_nonfinite_events.exchange(0)
    };
}

AMRRestrictionPPCounters
HermiteWENO2D::take_restriction_pp_counters () noexcept
{
    return {
        m_restrict_rho_events.exchange(0),
        m_restrict_pressure_events.exchange(0),
        m_restrict_any_events.exchange(0)
    };
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ncomp % sd_space_hweno == 0,
        "HermiteWENO2D interpolation requires complete SD-packed variables.");

    const Box target_fine_region = fine_region & fine.box();
    if (!target_fine_region.ok()) {
        return;
    }
    if (m_prolongation_pp_enabled) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ncomp / sd_space_hweno == NEQ,
            "HermiteWENO2D prolongation PP requires a complete SD-packed state.");
    }
    const Box parent_region = amrex::coarsen(target_fine_region, ratio);
    const Box full_fine_region = amrex::refine(parent_region, ratio);

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    const int nvar = ncomp / sd_space_hweno;
    const auto prolong_weight_mode = m_prolong_weight_mode;

    Box required_coarse_stencil = amrex::coarsen(full_fine_region, ratio);
    required_coarse_stencil.grow(IntVect(AMREX_D_DECL(1,1,0)));
    if (!crse.box().contains(required_coarse_stencil)) {
        const char* mode_name =
            prolong_weight_mode == ProlongWeightMode::Conservative
            ? "conservative" : "childwise";
        std::ostringstream message;
        message << "[HWENO_PROLONG_BOX_MISMATCH] rank="
                << ParallelDescriptor::MyProc()
                << " fine_region=" << fine_region
                << " target_fine_region=" << target_fine_region
                << " full_fine_region=" << full_fine_region
                << " crse_box=" << crse.box()
                << " required_coarse_stencil=" << required_coarse_stencil
                << " ratio=" << ratio
                << " crse_comp=" << crse_comp
                << " ncomp=" << ncomp
                << " weight_mode=" << mode_name;
        amrex::Abort(message.str());
    }

    Array4<Real const> const& carr = crse.const_array(crse_comp);
    Array4<Real>       const& farr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 2)
    Box by = amrex::coarsen(full_fine_region, IntVect(AMREX_D_DECL(ratio[0],1,1)));
    by.grow(IntVect(AMREX_D_DECL(1,0,0))); // halo for x-direction stage
    FArrayBox tmpy(by, ncomp);
#ifdef AMREX_USE_GPU
    Elixir tmpy_eli;
    if (run_on_gpu) { tmpy_eli = tmpy.elixir(); }
#endif
    Array4<Real> const& tmpyarr = tmpy.array();
    Array4<Real const> srcarr = carr;

    const Real hy = crse_geom.CellSize(1);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_interp_y(
            i, j, k, n, tmpyarr, srcarr, ratio, prolong_weight_mode, hy);
    });

    const Box intermediate_parent_region = amrex::coarsen(
        by, IntVect(AMREX_D_DECL(1,ratio[1],1)));
    const IntVect intermediate_child_ratio(AMREX_D_DECL(1,ratio[1],1));
    IArrayBox intermediate_nonfinite(intermediate_parent_region, 1);
    intermediate_nonfinite.setVal(0);
    if (m_prolongation_pp_enabled) {
        const auto y_counts = apply_parentwise_prolongation_pp(
            tmpy, crse, intermediate_parent_region, crse_comp,
            intermediate_child_ratio, nvar, true, false,
            m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, nullptr, runon);
        m_pp_rho_events.fetch_add(y_counts.rho);
        m_pp_pressure_events.fetch_add(y_counts.pressure);
        m_pp_any_events.fetch_add(y_counts.any);
        m_pp_nonfinite_events.fetch_add(y_counts.nonfinite);
    } else {
        const Long intermediate_nonfinite_count = hweno_check_intermediate_finite(
            tmpy, crse, intermediate_parent_region, intermediate_child_ratio,
            crse_comp, nvar, intermediate_nonfinite, runon);
        m_pp_nonfinite_events.fetch_add(intermediate_nonfinite_count);
    }

    srcarr = tmpy.const_array();
    FArrayBox raw_fine(full_fine_region, ncomp);
#ifdef AMREX_USE_GPU
    Elixir raw_fine_eli;
    if (run_on_gpu) { raw_fine_eli = raw_fine.elixir(); }
#endif
    auto const& raw = raw_fine.array();
    const Real hx = crse_geom.CellSize(0);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, full_fine_region, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_interp_x(
            i, j, k, n, raw, srcarr, ratio, prolong_weight_mode, hx);
    });

    const auto counts = apply_parentwise_prolongation_pp(
        raw_fine, crse, parent_region, crse_comp, ratio, nvar,
        m_prolongation_pp_enabled, true,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p,
        &intermediate_nonfinite, runon);
    m_pp_rho_events.fetch_add(counts.rho);
    m_pp_pressure_events.fetch_add(counts.pressure);
    m_pp_any_events.fetch_add(counts.any);
    m_pp_nonfinite_events.fetch_add(counts.nonfinite);

    auto const& raw_const = raw_fine.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp, i, j, k, n,
    {
        farr(i,j,k,n) = raw_const(i,j,k,n);
    });
#else
    amrex::ignore_unused(
        full_fine_region, run_on_gpu, nvar, carr, farr, crse_geom);
    amrex::Abort("HermiteWENO2D prolongation is implemented only in 2D.");
#endif
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
    if (m_prolongation_pp_enabled) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ncomp / sd_space_hweno == NEQ,
            "HermiteWENO2D restriction PP requires a complete SD-packed state.");
    }

    const Box target_crse_region = crse_region & crse.box();
    if (!target_crse_region.ok()) {
        return;
    }

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);
    Array4<Real const> const& finearr = fine.const_array(fine_comp);
    Array4<Real>       const& crsearr = crse.array(crse_comp);

#if (AMREX_SPACEDIM == 2)
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, target_crse_region, ncomp/sd_space_hweno, i, j, k, n,
    {
        Mortar2D::mortar_restrict_with_type(
            i, j, k, n, crsearr, finearr, ratio,
            Mortar2D::Type::ScaleRef);
    });
#else
    amrex::ignore_unused(fine_geom, crse_geom);
#endif

    const Box inner_box = amrex::grow(target_crse_region, -1);
    if (!inner_box.ok()) {
        const auto counts = apply_restriction_positivity(
            crse, target_crse_region, crse_comp, ncomp/sd_space_hweno,
            &fine, fine_comp, ratio,
            m_prolongation_pp_enabled,
            m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, runon);
        m_restrict_rho_events.fetch_add(counts.rho);
        m_restrict_pressure_events.fetch_add(counts.pressure);
        m_restrict_any_events.fetch_add(counts.any);
        return;
    }

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::refine(inner_box, IntVect(ratio[0],ratio[1],1));
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
    Box by = amrex::refine(
        inner_box, IntVect(AMREX_D_DECL(ratio[0],1,1)));
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
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_restrict_y(i, j, k, n, tmpyarr, srcarr, ratio, hy_f, hy_c);
    });
    const IntVect y_source_ratio(AMREX_D_DECL(1,ratio[1],1));
#if (AMREX_SPACEDIM == 2)
    FArrayBox const* y_source = &fine;
    const int y_source_comp = fine_comp;
#else
    FArrayBox const* y_source = &tmpz;
    const int y_source_comp = 0;
#endif
    const auto y_counts = apply_restriction_positivity(
        tmpy, by, 0, ncomp/sd_space_hweno,
        y_source, y_source_comp, y_source_ratio,
        m_prolongation_pp_enabled,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, runon);
    m_restrict_rho_events.fetch_add(y_counts.rho);
    m_restrict_pressure_events.fetch_add(y_counts.pressure);
    m_restrict_any_events.fetch_add(y_counts.any);
    amr_assert_fab_finite(
        tmpy, by, 0, ncomp, runon,
        "HWENO restriction produced a non-finite intermediate state.");
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = finearr;
#else
    srcarr = tmpy.const_array();
#endif
    const Real hx_f = fine_geom.CellSize(0);
    const Real hx_c = crse_geom.CellSize(0);
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, inner_box, ncomp/sd_space_hweno, i, j, k, n,
    {
        hweno_restrict_x(i, j, k, n, crsearr, srcarr, ratio, hx_f, hx_c);
    });
    const IntVect x_source_ratio(AMREX_D_DECL(ratio[0],1,1));
    const auto x_counts = apply_restriction_positivity(
        crse, inner_box, crse_comp, ncomp/sd_space_hweno,
        &tmpy, 0, x_source_ratio,
        m_prolongation_pp_enabled,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, runon);
    m_restrict_rho_events.fetch_add(x_counts.rho);
    m_restrict_pressure_events.fetch_add(x_counts.pressure);
    m_restrict_any_events.fetch_add(x_counts.any);
    amr_assert_fab_finite(
        crse, target_crse_region, crse_comp, ncomp, runon,
        "HWENO restriction produced a non-finite coarse polynomial.");
    const auto counts = apply_restriction_positivity(
        crse, target_crse_region, crse_comp, ncomp/sd_space_hweno,
        &fine, fine_comp, ratio,
        m_prolongation_pp_enabled,
        m_pp_gamma, m_pp_eps_rho, m_pp_eps_p, runon);
    m_restrict_rho_events.fetch_add(counts.rho);
    m_restrict_pressure_events.fetch_add(counts.pressure);
    m_restrict_any_events.fetch_add(counts.any);
}

}
