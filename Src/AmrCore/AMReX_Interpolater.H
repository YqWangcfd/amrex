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
        {  0.603383141696693,   0.771704926696188,  -0.856836867851174,   2.469865488345720, -12.313391031781396,  11.580862183718867, -1.281980197420377,   0.073293286821232,  -0.107075962763792,   0.060175032538039},
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

    // ## 2D mortar : scale refinement
    // ## fine -> coarse projection matrices
    static constexpr Real Px2D[2][5][5] = {
    // ---------- Px2D[0] ----------
    {
        {  0.641042236219882,         0.618076982602313,        -0.332582031664555,       -0.0466770180530059,        0.0166769784883341},
        {-0.0744899769895421,         0.380675785616308,         0.652619642561019,         0.151492669414719,      -0.00583818316292594},
        {  0.042709608237357,        -0.104991111752765,                    0.0375,         0.378167386330886,         0.146614117184522},
        {-0.0399601549200435,        0.0767502392044148,       0.00762425534322947,        -0.126429086614335,       -0.0224451904528431},
        { 0.0613301128925325,        -0.100919064504609,       -0.0448493662396924,         0.151940105480656,        0.0359610647781447}
    },
    // ---------- Px2D[1] ----------
    {
        {  0.0359610647781447,         0.151940105480656,       -0.0448493662396923,        -0.100919064504609,        0.0613301128925323},
        { -0.0224451904528431,        -0.126429086614335,       0.00762425534322946,        0.0767502392044148,       -0.0399601549200435},
        {   0.146614117184522,         0.378167386330886,                    0.0375,        -0.104991111752765,         0.042709608237357},
        {-0.00583818316292593,         0.151492669414719,         0.652619642561019,         0.380675785616308,        -0.074489976989542},
        {  0.0166769784883341,       -0.0466770180530058,        -0.332582031664556,         0.618076982602314,         0.641042236219883}
    }
    };

    static constexpr Real Py2D[2][5][5] = {
    // ---------- Py2D[0] ----------
    {
        {  0.641042236219882,       -0.0744899769895421,         0.042709608237357,       -0.0399601549200435,        0.0613301128925325},
        {  0.618076982602313,         0.380675785616308,        -0.104991111752765,        0.0767502392044148,        -0.100919064504609},
        { -0.332582031664555,         0.652619642561019,                    0.0375,       0.00762425534322947,       -0.0448493662396924},
        {-0.0466770180530059,         0.151492669414719,         0.378167386330886,        -0.126429086614335,         0.151940105480656},
        { 0.0166769784883341,      -0.00583818316292594,         0.146614117184522,       -0.0224451904528431,        0.0359610647781447}
    },
    // ---------- Py2D[1] ----------
    {
        { 0.0359610647781447,       -0.0224451904528431,         0.146614117184522,      -0.00583818316292593,        0.0166769784883341},
        {  0.151940105480656,        -0.126429086614335,         0.378167386330886,         0.151492669414719,       -0.0466770180530058},
        {-0.0448493662396923,       0.00762425534322946,                    0.0375,         0.652619642561019,        -0.332582031664556},
        { -0.100919064504609,        0.0767502392044148,        -0.104991111752765,         0.380675785616308,         0.618076982602314},
        { 0.0613301128925323,       -0.0399601549200435,         0.042709608237357,        -0.074489976989542,         0.641042236219883}
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

    // ## 2D mortar : scale refinement
    // ## coarse -> fine projection matrices
    static constexpr Real Px2D_prol[2][5][5] = {
    // ---------- Px2D[0] ----------
    {
        {  0.774190504909649,         0.315598720481656,        -0.135925017591863,        0.0665826964702072,       -0.0204469042696481},
        {-0.0615480890323645,         0.303137941289733,         0.879740570089535,        -0.163695536273256,        0.0423651139263522},
        { 0.0316768880649073,        -0.101905089898884,                       0.2,        -0.392522101101031,          1.26275030293501},
        {   11.5584026375638,          -34.238314463223,           53.754938531917,         -59.5878083897827,          29.5127816835248},
        {   45.2224569037659,         -130.717816290963,          194.301245915585,         -197.116366471972,          89.3104799435833}
    },
    // ---------- Px2D[1] ----------
    {
        {   89.3104799435833,         -197.116366471972,          194.301245915585,         -130.717816290963,          45.2224569037659},
        {   29.5127816835248,         -59.5878083897827,           53.754938531917,          -34.238314463223,          11.5584026375638},
        {   1.26275030293501,         -0.39252210110103,         0.200000000000001,        -0.101905089898886,        0.0316768880649077},
        { 0.0423651139263526,        -0.163695536273256,         0.879740570089533,         0.303137941289734,        -0.061548089032365},
        {-0.0204469042696493,        0.0665826964702087,         -0.13592501759186,         0.315598720481653,         0.774190504909651}
    }
    };

    static constexpr Real Py2D_prol[2][5][5] = {
    // ---------- Py2D[0] ----------
    {
        {  0.774190504909649,       -0.0615480890323645,        0.0316768880649073,          11.5584026375638,          45.2224569037659},
        {  0.315598720481656,         0.303137941289733,        -0.101905089898884,          -34.238314463223,         -130.717816290963},
        { -0.135925017591863,         0.879740570089535,                       0.2,           53.754938531917,          194.301245915585},
        { 0.0665826964702072,        -0.163695536273256,        -0.392522101101031,         -59.5878083897827,         -197.116366471972},
        {-0.0204469042696481,        0.0423651139263522,          1.26275030293501,          29.5127816835248,          89.3104799435833}
    },
    // ---------- Py2D[1] ----------
    {
        { 89.3104799435833,          29.5127816835248,          1.26275030293501,        0.0423651139263526,       -0.0204469042696493},
        {-197.116366471972,         -59.5878083897827,         -0.39252210110103,        -0.163695536273256,        0.0665826964702087},
        { 194.301245915585,           53.754938531917,         0.200000000000001,         0.879740570089533,         -0.13592501759186},
        {-130.717816290963,          -34.238314463223,        -0.101905089898886,         0.303137941289734,         0.315598720481653},
        { 45.2224569037659,          11.5584026375638,        0.0316768880649077,        -0.061548089032365,         0.774190504909651}
    }
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
