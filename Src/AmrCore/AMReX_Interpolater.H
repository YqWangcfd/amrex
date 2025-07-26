#ifndef AMREX_INTERPOLATER_H_
#define AMREX_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_InterpBase.H>
#include <SD_Mapping_K.H>

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
* \brief 2D Mortar Prolongation and Restriction
*/
class Mortar2D
    :
    public Interpolater
{
public:
    static constexpr int sd_order  = 5;
    static constexpr int sd_space  = 25;

    enum class Type : int { OrderRef, ScaleRef };
    Type type = Type::OrderRef;

    // ------------------------------------------------------------------------------------------------
    // ## 2D mortar : order refinement
    // ## fine -> coarse projection matrices
    // ----------------------------------------------------------------------------------------------
    static constexpr Real Px[5][10] = {
        { 0.603383141696693,   0.771704926696188,  -0.856836867851174,   2.469865488345720, -12.313391031781396,  11.580862183718867, -1.281980197420377,   0.073293286821232,  -0.107075962763792,   0.060175032538039},
        {-0.0637108487698427,   0.338463378171787,   0.783859428088837,  -0.390843753679491,   2.266194771743281,  -2.026752144784240,  0.044770809223417,   0.018834799716061,   0.065614679953372,  -0.036431119663181},
        { 0.0384207195631531,  -0.0887505790994087,  -0.006983256154527,   0.496287373858870,   0.061025741831907,   0.061025741831921,  0.496287373858867,  -0.006983256154527,  -0.088750579099409,   0.038420719563153},
        {-0.0364311196631815,   0.0656146799533715,   0.018834799716061,   0.0447708092234142,  -2.026752144784220,   2.266194771743262, -0.390843753679487,   0.783859428088836,   0.338463378171787,  -0.063710848769843},
        { 0.0601750325380391,  -0.107075962763791,   0.073293286821232,  -1.281980197420372,  11.580862183718832, -12.313391031781363,  2.469865488345713,  -0.856836867851173,   0.771704926696187,   0.603383141696693}
    };

    static constexpr Real Py[10][5] = {
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

    static constexpr Real Px2D[2][5][5] = {
    // ---------- Px2D[0] ----------
    {
        {  44.6552399717917,  -98.5581832359859,   97.1506229577926,  -65.3589081454813,   22.611228451883  },
        {  14.7563908417624,  -29.7939041948913,   26.8774692659585,  -17.1191572316115,    5.77920131878192},
        {   0.631375151467504,-0.196261050550514,  0.0999999999999995,-0.0509525449494418,  0.0158384440324533},
        {  0.0211825569631767,-0.0818477681366283,  0.439870285044767,  0.151568970644867, -0.0307740445161821},
        { -0.0102234521348255, 0.0332913482351046, -0.0679625087959305, 0.157799360240828,  0.387095252454824}
    },
    // ---------- Px2D[1] ----------
    {
        {   0.387095252454824,  0.157799360240827, -0.0679625087959325, 0.0332913482351036, -0.0102234521348257},
        { -0.0307740445161820,  0.151568970644867,  0.439870285044768,-0.0818477681366281,  0.0211825569631771},
        {  0.0158384440324533,-0.0509525449494419,  0.0999999999999982,-0.196261050550515,   0.631375151467503},
        {   5.77920131878192, -17.1191572316115,    26.8774692659585, -29.7939041948913,    14.7563908417624},
        {  22.6112284518829,  -65.3589081454813,    97.1506229577926, -98.5581832359859,    44.6552399717917}
    }
    };

    static constexpr Real Py2D[2][5][5] = {
    // ---------- Py2D[0] ----------
    {
        { 44.6552399717917,          14.7563908417624,         0.631375151467504,        0.0211825569631767,       -0.0102234521348255},
        {-98.5581832359859,         -29.7939041948913,        -0.196261050550514,       -0.0818477681366283,        0.0332913482351046},
        { 97.1506229577926,          26.8774692659585,        0.0999999999999995,         0.439870285044767,       -0.0679625087959305},
        {-65.3589081454813,         -17.1191572316115,       -0.0509525449494418,         0.151568970644867,         0.157799360240828},
        {  22.611228451883,          5.77920131878192,        0.0158384440324533,       -0.0307740445161821,         0.387095252454824}
    },
    // ---------- Py2D[1] ----------
    {
        {  0.387095252454824,        -0.030774044516182,        0.0158384440324533,          5.77920131878192,          22.6112284518829},
        {  0.157799360240827,         0.151568970644867,       -0.0509525449494419,         -17.1191572316115,         -65.3589081454813},
        {-0.0679625087959325,         0.439870285044768,        0.0999999999999982,          26.8774692659585,          97.1506229577926},
        { 0.0332913482351036,       -0.0818477681366281,        -0.196261050550515,         -29.7939041948913,         -98.5581832359859},
        {-0.0102234521348257,        0.0211825569631771,         0.631375151467503,          14.7563908417624,          44.6552399717917}
    }
    };


    // ------------------------------------------------------------------------------------------------
    // ## 2D mortar : order refinement
    // ## coarse -> fine projection matrices
    // ----------------------------------------------------------------------------------------------
    static constexpr Real Px_prol[10][5] = {
        {    1.12654974652665,        -0.186142938148128,        0.0914519644157589,       -0.0461683609041033,        0.0143095881098134},
        {   0.390101138933886,         0.778774529455574,        -0.249911278368893,         0.116171197060536,       -0.0351355870811183},
        { -0.0685095954079369,         0.921586454726534,          0.20000000000003,        -0.074372859226575,        0.0212959999079791},
        { -0.0714627134813877,         0.365161025548818,         0.835853630900096,        -0.175557601783677,        0.0460056588160938},
        {-0.00814181071916837,        0.0350615631541298,         0.997605683052834,       -0.0322587612893716,       0.00773332580144883},
        { 0.00773332580144507,       -0.0322587612893602,         0.997605683052832,        0.0350615631541189,      -0.00814181071916481},
        {  0.0460056588160917,        -0.175557601783671,         0.835853630900083,         0.365161025548815,       -0.0714627134813867},
        {  0.0212959999079802,       -0.0743728592265789,         0.200000000000025,         0.921586454726536,       -0.0685095954079378},
        { -0.0351355870811189,         0.116171197060538,        -0.249911278368887,         0.778774529455573,         0.390101138933886},
        {  0.0143095881098131,       -0.0461683609041024,        0.0914519644157527,        -0.186142938148129,          1.12654974652665}
    };

    static constexpr Real Py_prol[5][10] = {
        {   1.12654974652665,         0.390101138933886,       -0.0685095954079369,       -0.0714627134813877,      -0.00814181071916837,       0.00773332580144507,        0.0460056588160917,        0.0212959999079802,       -0.0351355870811189,        0.0143095881098131},
        { -0.186142938148128,         0.778774529455574,         0.921586454726534,         0.365161025548818,        0.0350615631541298,       -0.0322587612893602,        -0.175557601783671,       -0.0743728592265789,         0.116171197060538,       -0.0461683609041024},
        { 0.0914519644157589,        -0.249911278368893,          0.20000000000003,         0.835853630900096,         0.997605683052834,         0.997605683052832,         0.835853630900083,         0.200000000000025,        -0.249911278368887,        0.0914519644157527},
        {-0.0461683609041033,         0.116171197060536,        -0.074372859226575,        -0.175557601783677,       -0.0322587612893716,        0.0350615631541189,         0.365161025548815,         0.921586454726536,         0.778774529455573,        -0.186142938148129},
        { 0.0143095881098134,       -0.0351355870811183,        0.0212959999079791,        0.0460056588160938,       0.00773332580144883,      -0.00814181071916481,       -0.0714627134813867,       -0.0685095954079378,         0.390101138933886,          1.12654974652665}
    };

    Mortar2D() = default;              
    explicit Mortar2D(Type t) noexcept : type(t) {}

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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mortar_interp(const int i, const int j, const int k, const int n, 
                Array4<Real> const& finearr, Array4<const Real> const& crsearr,
                const IntVect&   ratio);

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


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mortar_restrict(const int i, const int j, const int k, const int n, 
                Array4<Real> const crsearr, Array4<const Real> const& finearr,
                const IntVect&   ratio);
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
extern AMREX_EXPORT Mortar2D                  mortar_interp_orderRef;
extern AMREX_EXPORT Mortar2D                  mortar_interp_scaleRef;

}

#endif /*_INTERPOLATER_H_*/
