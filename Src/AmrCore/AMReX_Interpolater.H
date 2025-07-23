#ifndef AMREX_INTERPOLATER_H_
#define AMREX_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_InterpBase.H>

namespace amrex {

class Geometry;
class FArrayBox;
class IArrayBox;

/**
* \brief Virtual base class for interpolaters.
*
* Specifies interpolater interface for coarse-to-fine interpolation in space.
*/
class Interpolater
    : public InterpBase
{
public:
    /**
    * \brief Coarse to fine interpolation in space.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const & bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) = 0;

    virtual void Lagrange_interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const & bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon)
    {
        amrex::Abort("Interpolater::Lagrange_interp not implemented");
    }

    virtual void restrict (const FArrayBox& fine,
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
        amrex::Abort("Interpolater::restrict not implemented");
    }

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param fine_values
    * \param fine_known
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    virtual void interp_face (const FArrayBox& /*crse*/,
                              const int        /*crse_comp*/,
                              FArrayBox&       /*fine*/,
                              const int        /*fine_comp*/,
                              const int        /*ncomp*/,
                              const Box&       /*fine_region*/,
                              const IntVect&   /*ratio*/,
                              const IArrayBox& /*solve_mask*/,
                              const Geometry&  /*crse_geom*/,
                              const Geometry&  /*fine_geom*/,
                              Vector<BCRec> const & /*bcr*/,
                              const int        /*bccomp*/,
                              RunOn            /*runon*/)
    { amrex::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead."); }

    /**
    * \brief Coarse to fine interpolation in space.
    *        Takes an Array of FArrayBox*s
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param fine_values
    * \param fine_known
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    virtual void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& /*crse*/,
                             const int         /*crse_comp*/,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& /*fine*/,
                             const int         /*fine_comp*/,
                             const int         /*ncomp*/,
                             const Box&        /*fine_region*/,
                             const IntVect&    /*ratio*/,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& /*solve_mask*/,
                             const Geometry&   /*crse_geom*/,
                             const Geometry&   /*fine_geom*/,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                             const int         /*actual_comp*/,
                             const int         /*actual_state*/,
                             const RunOn       /*runon*/)

    { amrex::Abort("The Array<FArrayBox*, AMREX_SPACEDIM> version of this Interpolater is not implemented or does not apply. Call 'interp' instead."); }

    /**
    * \brief Re-visit the interpolation to protect against under- or overshoots.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param fine_state
    * \param state_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    */
    virtual void protect (const FArrayBox& /*crse*/,
                          int              /*crse_comp*/,
                          FArrayBox&       /*fine*/,
                          int              /*fine_comp*/,
                          FArrayBox&       /*fine_state*/,
                          int              /*state_comp*/,
                          int              /*ncomp*/,
                          const Box&       /*fine_region*/,
                          const IntVect&   /*ratio*/,
                          const Geometry&  /*crse_geom*/,
                          const Geometry&  /*fine_geom*/,
                          Vector<BCRec>&   /*bcr*/,
                          RunOn            /*runon*/) {}
};


/**
* \brief Bilinear interpolation on node centered data.
*
* Bilinear interpolation on node centered data.
*/
class NodeBilinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;
};


/**
* \brief Bilinear interpolation on cell centered data.
*
* Bilinear interpolation on cell centered data.
*/
class CellBilinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;
};


/**
* \brief Linear conservative interpolation on cell centered data
*
*
* Linear conservative interpolation on cell centered data
* I.e, conservative interpolation with a limiting scheme that
* preserves the value of any linear combination  of the
* fab components; //!< e.g.,
*
* if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
* sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
* in all fine cells if,jf covering coarse cell ic,jc.
*/
class CellConservativeLinear
    :
    public Interpolater
{
public:

    /**
    * \brief The constructor.
    *
    * \param do_linear_limiting_
    */
    explicit CellConservativeLinear (bool do_linear_limiting_= true);

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;
    /**
    * \brief Coarse to fine interpolation in space.
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

protected:

    bool do_linear_limiting;
};


/**
* \brief Lin. cons. interp. on cc data with protection against under/over-shoots.
*
* Linear conservative interpolation on cell centered data
* but with protection against undershoots or overshoots.
*/
class CellConservativeProtected
    :
    public CellConservativeLinear
{
public:

    /**
    * \brief The constructor.
    */
    CellConservativeProtected ();

    /**
    * \brief Re-visit the interpolation to protect against under- or overshoots.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param fine_state
    * \param state_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    */
    void protect (const FArrayBox& crse,
                  int              crse_comp,
                  FArrayBox&       fine,
                  int              fine_comp,
                  FArrayBox&       fine_state,
                  int              state_comp,
                  int              ncomp,
                  const Box&       fine_region,
                  const IntVect&   ratio,
                  const Geometry&  crse_geom,
                  const Geometry&  fine_geom,
                  Vector<BCRec>&   bcr,
                  RunOn            runon) override;
};


/**
* \brief Quadratic interpolation on cell centered data.
*
* Quadratic interpolation on cell centered data.
*/
class CellQuadratic
    :
    public Interpolater
{
public:

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;
};


/**
* \brief Piecewise Constant interpolation on cell centered data.
*/
class PCInterp
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                 int              crse_comp,
                 FArrayBox&       fine,
                 int              fine_comp,
                 int              ncomp,
                 const Box&       fine_region,
                 const IntVect&   ratio,
                 const Geometry&  crse_geom,
                 const Geometry&  fine_geom,
                 Vector<BCRec> const&  bcr,
                 int              actual_comp,
                 int              actual_state,
                 RunOn            runon) override;
};


/**
* \brief Conservative quartic interpolation on cell averaged data.
*
* An order 4 polynomial is used to fit the data.  For each cell involved
* in constructing the polynomial, the average of the polynomial inside that
* cell is equal to the cell averaged value of the original data.
*/
class CellConservativeQuartic
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
                 int              crse_comp,
                 FArrayBox&       fine,
                 int              fine_comp,
                 int              ncomp,
                 const Box&       fine_region,
                 const IntVect&   ratio,
                 const Geometry&  crse_geom,
                 const Geometry&  fine_geom,
                 Vector<BCRec> const&  bcr,
                 int              actual_comp,
                 int              actual_state,
                 RunOn            runon) override;
};

/**
* \brief Divergence-preserving interpolation on face centered data.
*
* Divergence-preserving interpolation on face centered data,
* i.e., it ensures the divergence of the fine ghost cells match the value
* of the divergence of the underlying crse cell. All fine cells overlying
* a given coarse cell will have the same divergence, even when the coarse
* grid divergence is spatially varying.
* Based on Vanella et. al. (doi:10.1016/j.jcp.2010.05.003, section 3.2),
* but solves the interior closure problem using least squares with an
* initial guess equal to the average of fine face values across the cell.
*/
class FaceDivFree
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    * \param runon
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param solve_mask
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    * \param runon
    */
    void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                     int         crse_comp,
                     Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                     int         fine_comp,
                     int         ncomp,
                     const Box&        fine_region,
                     const IntVect&    ratio,
                     Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                     const Geometry&   crse_geom,
                     const Geometry&   fine_geom,
                     Vector<Array<BCRec, AMREX_SPACEDIM> > const& bcr,
                     int         actual_comp,
                     int         actual_state,
                     RunOn       runon) override;
};


/**
* \brief Piecewise constant tangential interpolation / linear normal interpolation of face data.
*
* Piecewise constant tangential interpolation / linear normal interpolation of face data.
*/
class FaceLinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param solve_mask
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    void interp_face (const FArrayBox& crse,
                      int        crse_comp,
                      FArrayBox&       fine,
                      int        fine_comp,
                      int        ncomp,
                      const Box&       fine_region,
                      const IntVect&   ratio,
                      const IArrayBox& solve_mask,
                      const Geometry&  crse_geom,
                      const Geometry&  fine_geom,
                      Vector<BCRec> const & bcr,
                      int        bccomp,
                      RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                     int         crse_comp,
                     Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                     int         fine_comp,
                     int         ncomp,
                     const Box&        fine_region,
                     const IntVect&    ratio,
                     Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                     const Geometry&   /*crse_geom*/,
                     const Geometry&   /*fine_geom*/,
                     Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                     int         /*actual_comp*/,
                     int         /*actual_state*/,
                     RunOn       runon) override;


};

/**
* \brief Bilinear tangential interpolation / linear normal interpolation of face data.
*
* Bilinear tangential interpolation / linear normal interpolation of face data.
*/
class FaceConservativeLinear
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space for face-based data.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param solve_mask
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param bccomp
    * \param runon
    */
    void interp_face (const FArrayBox& crse,
                      int        crse_comp,
                      FArrayBox&       fine,
                      int        fine_comp,
                      int        ncomp,
                      const Box&       fine_region,
                      const IntVect&   ratio,
                      const IArrayBox& solve_mask,
                      const Geometry&  crse_geom,
                      const Geometry&  fine_geom,
                      Vector<BCRec> const & bcr,
                      int        bccomp,
                      RunOn            runon) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                     int         crse_comp,
                     Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                     int         fine_comp,
                     int         ncomp,
                     const Box&        fine_region,
                     const IntVect&    ratio,
                     Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                     const Geometry&   /*crse_geom*/,
                     const Geometry&   /*fine_geom*/,
                     Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                     int         /*actual_comp*/,
                     int         /*actual_state*/,
                     RunOn       runon) override;


};

/**
* \brief Quartic interpolation on cell centered data.
*
* Quartic interpolation on cell centered data.
*/
class CellQuartic
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    // void Lagrange_interp (const FArrayBox& crse,
    //              int              crse_comp,
    //              FArrayBox&       fine,
    //              int              fine_comp,
    //              int              ncomp,
    //              const Box&       fine_region,
    //              const IntVect&   ratio,
    //              const Geometry&  crse_geom,
    //              const Geometry&  fine_geom,
    //              Vector<BCRec> const& bcr,
    //              int              actual_comp,
    //              int              actual_state,
    //              RunOn            runon) override;
    /**
    * \brief Fine to coarse interpolation in space.
    *
    * \param fine
    * \param fine_comp
    * \param crse
    * \param crse_comp
    * \param ncomp
    * \param crse_region
    * \param ratio
    * \param fine_geom
    * \param crse_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void restrict (const FArrayBox& fine,
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
                 RunOn            runon) override;
};


/**
* \brief Central WENO interpolation on cell centered data.
*
* Central WENO interpolation on cell centered data.
*/
class CellWENO
    :
    public Interpolater
{
public:
    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    /**
    * \brief Fine to coarse interpolation in space.
    *
    * \param fine
    * \param fine_comp
    * \param crse
    * \param crse_comp
    * \param ncomp
    * \param crse_region
    * \param ratio
    * \param fine_geom
    * \param crse_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void restrict (const FArrayBox& fine,
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
                 RunOn            runon) override;
};


/**
* \brief Central WENO interpolation on cell centered data.
*
* Central WENO interpolation on cell centered data.
*/
class 2DMortar
    :
    public Interpolater
{
public:
    constexpr int sd_order  = 5;
    constexpr int sd_space  = 25;
    // ------------------------------------------------------------------------------------------------
    // ## 2D mortar : order refinement
    // ## fine -> coarse projection matrices
    // ----------------------------------------------------------------------------------------------
    Real Px = [5][10] = {
        { 0.603383141696693,   0.771704926696189,  -0.856836867851174,   2.469865488345720, -12.313391031781396,  11.580862183718867, -1.281980197420377,   0.073293286821232,  -0.107075962763792,   0.060175032538039},
        {-0.063710848769843,   0.338463378171787,   0.783859428088837,  -0.390843753679491,   2.266194771743281,  -2.026752144784240,  0.044770809223417,   0.018834799716061,   0.065614679953372,  -0.036431119663181},
        { 0.038420719563153,  -0.088750579099409,  -0.006983256154527,   0.496287373858870,   0.061025741831907,   0.061025741831921,  0.496287373858867,  -0.006983256154527,  -0.088750579099409,   0.038420719563153},
        {-0.036431119663181,   0.065614679953372,   0.018834799716061,   0.044770809223414,  -2.026752144784220,   2.266194771743262, -0.390843753679487,   0.783859428088836,   0.338463378171787,  -0.063710848769843},
        { 0.060175032538039,  -0.107075962763791,   0.073293286821232,  -1.281980197420372,  11.580862183718832, -12.313391031781363,  2.469865488345713,  -0.856836867851173,   0.771704926696187,   0.603383141696693}
    };

    Real Py = [10][5] = {
        { 0.603383141696693,  -0.063710848769843,   0.038420719563153,  -0.036431119663181,   0.060175032538039},
        { 0.771704926696189,   0.338463378171787,  -0.088750579099409,   0.065614679953372,  -0.107075962763791},
        {-0.856836867851174,   0.783859428088837,  -0.006983256154527,   0.018834799716061,   0.073293286821232},
        { 2.469865488345720,  -0.390843753679491,   0.496287373858870,   0.044770809223414,  -1.281980197420372},
        {-12.313391031781396,   2.266194771743281,   0.061025741831907,  -2.026752144784220,  11.580862183718832},
        {11.580862183718867,  -2.026752144784240,   0.061025741831921,   2.266194771743262, -12.313391031781363},
        {-1.281980197420377,   0.044770809223417,   0.496287373858867,  -0.390843753679487,   2.469865488345713},
        { 0.073293286821232,   0.018834799716061,  -0.006983256154527,   0.783859428088836,  -0.856836867851173},
        {-0.107075962763792,   0.065614679953372,  -0.088750579099409,   0.338463378171787,   0.771704926696187},
        { 0.060175032538039,  -0.036431119663181,   0.038420719563153,  -0.063710848769843,   0.603383141696693}
    }; 

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, int ratio) override;

    /**
    * \brief Returns coarsened box given fine box and refinement ratio.
    *
    * \param fine
    * \param ratio
    */
    Box CoarseBox (const Box& fine, const IntVect& ratio) override;

    /**
    * \brief Coarse to fine interpolation in space.
    *
    * \param crse
    * \param crse_comp
    * \param fine
    * \param fine_comp
    * \param ncomp
    * \param fine_region
    * \param ratio
    * \param crse_geom
    * \param fine_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void interp (const FArrayBox& crse,
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
                 RunOn            runon) override;

    /**
    * \brief Fine to coarse interpolation in space.
    *
    * \param fine
    * \param fine_comp
    * \param crse
    * \param crse_comp
    * \param ncomp
    * \param crse_region
    * \param ratio
    * \param fine_geom
    * \param crse_geom
    * \param bcr
    * \param actual_comp
    * \param actual_state
    */
    void restrict (const FArrayBox& fine,
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
                 RunOn            runon) override;
};

//! CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
extern AMREX_EXPORT PCInterp                  pc_interp;
extern AMREX_EXPORT NodeBilinear              node_bilinear_interp;
extern AMREX_EXPORT FaceDivFree               face_divfree_interp;
extern AMREX_EXPORT FaceLinear                face_linear_interp;
extern AMREX_EXPORT FaceConservativeLinear    face_cons_linear_interp;
extern AMREX_EXPORT CellConservativeLinear    lincc_interp;
extern AMREX_EXPORT CellConservativeLinear    cell_cons_interp;
extern AMREX_EXPORT CellBilinear              cell_bilinear_interp;
extern AMREX_EXPORT CellConservativeProtected protected_interp;
extern AMREX_EXPORT CellConservativeQuartic   quartic_interp;
extern AMREX_EXPORT CellQuadratic             quadratic_interp;
extern AMREX_EXPORT CellQuartic               cell_quartic_interp;
extern AMREX_EXPORT CellWENO                  cell_weno_interp;
extern AMREX_EXPORT 2DMortar                  mortar_interp_2d;

}

#endif /*_INTERPOLATER_H_*/
