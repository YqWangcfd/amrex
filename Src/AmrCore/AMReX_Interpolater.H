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
    // for 1D flattening of 2D L2-projection
    static constexpr Real P2D[25][100] = {
        {0.3640712156837722, 0.4656337431327621, -0.5170009212455949, 1.490275197926275, -7.429692565696124, 6.987697007968698, -0.7735252391124514, 0.04422393366747168, -0.06460783081261448, 0.03630860018450277, 0.4656337431327622, 0.5955284938871681, -0.6612252322956798, 1.906007365633272, -9.502304523562326, 8.937008402565388, -0.9893104342762522, 0.05656079053370133, -0.08263104799555519, 0.04643736907370814, -0.517000921245595, -0.6612252322956796, 0.7341694181090076, -2.116271809047847, 10.55056740429825, -9.92290968051373, 1.0984478970049, -0.06280039031442183, 0.09174663255667566, -0.05156018640273579, 1.490275197926275, 1.906007365633272, -2.116271809047848, 6.100235530521209, -30.41241955390237, 28.60317183285511, -3.166318646351201, 0.1810245596471841, -0.2644632250616785, 0.1486242361257831, -7.42969256569612, -9.502304523562326, 10.55056740429825, -30.41241955390239, 151.6195987015536, -142.5996845532993, 15.7855234658373, -0.9024889006343342, 1.318468199615021, -0.7409587059910421, 6.987697007968698, 8.937008402565388, -9.92290968051373, 28.60317183285511, -142.5996845532993, 134.1163689182889, -14.84643598858199, 0.8487994536684611, -1.240031967956475, 0.696878758723827, -0.7735252391124515, -0.9893104342762516, 1.0984478970049, -3.166318646351199, 15.7855234658373, -14.84643598858199, 1.643473226577976, -0.09396054230867072, 0.1372692638829015, -0.07714320009289266, 0.0442239336674717, 0.05656079053370128, -0.06280039031442186, 0.1810245596471841, -0.9024889006343335, 0.8487994536684609, -0.09396054230867078, 0.00537190589305935, -0.007847949250506081, 0.004410425919287459, -0.06460783081261448, -0.08263104799555521, 0.09174663255667564, -0.2644632250616786, 1.318468199615021, -1.240031967956475, 0.1372692638829016, -0.007847949250506081, 0.0114652618017928, -0.006443299543353001, 0.03630860018450277, 0.04643736907370813, -0.0515601864027358, 0.1486242361257831, -0.7409587059910417, 0.696878758723827, -0.07714320009289266, 0.004410425919287458, -0.006443299543352999, 0.003621034540954062},
        {-0.0384420520909106, 0.2042230964705688, 0.4729675643688154, -0.2358285320076577, 1.367383721071069, -1.222908076560414, 0.02701395152552696, 0.01136460062590503, 0.03959079173168834, -0.02198192343789861, -0.04916597587968335, 0.2611938564414028, 0.6049081824934112, -0.3016160502828915, 1.748833670207416, -1.564054615322047, 0.03454985404988384, 0.01453490773422041, 0.05063517178361035, -0.02811407452913547, 0.05458980410809174, -0.2900079008350409, -0.67163965719925, 0.3348893377219272, -1.941759230161198, 1.736595959647557, -0.03836127994615263, -0.01613835079531408, -0.05622107685630393, 0.03121552646451166, -0.1573572266098474, 0.8359590168154006, 1.936027349151025, -0.9653314985484587, 5.59719625659817, -5.005805175833199, 0.1105778765862207, 0.04651942179860329, 0.1620594335456815, -0.08997996515788478, 0.7844965938697592, -4.167631925366904, -9.651967652006341, 4.812611971384761, -27.90454237865326, 24.95619168322964, -0.5512804807771827, -0.2319202539091466, -0.8079392116910481, 0.448590622138372, -0.7378265592113007, 3.919697736843375, 9.07776800810548, -4.526307646729492, 26.24448933302283, -23.47153726930261, 0.5184845714699282, 0.218123219769652, 0.7598745657688141, -0.4219037760178737, 0.08167604648378236, -0.4339033483682339, -1.004892264371148, 0.5010539525025485, -2.90521682087244, 2.598256114692637, -0.05739529084690239, -0.02414584025836938, -0.08411672036029773, 0.04670397397805056, -0.004669577512512198, 0.02480709345482783, 0.05745163389044182, -0.02864622334071752, 0.1660968633981546, -0.1485473262632164, 0.003281399761630344, 0.001380464377809735, 0.004809115557505774, -0.002670156502692178, 0.006821900470529212, -0.03624129207802922, -0.0839325029340869, 0.04184997121544518, -0.2426549869946797, 0.2170164371863492, -0.0047938775013111, -0.002016754313060439, -0.007025755027445282, 0.003900897212498039, -0.003833802397751372, 0.02036704479442191, 0.04716876659049444, -0.02351903559495249, 0.1363683441271847, -0.121959876258931, 0.002694084901773292, 0.001133384685761435, 0.003948365501167152, -0.00219224381112914},
        {0.02318241447626296, -0.05355060324440214, -0.004213579037791132, 0.2994514348233649, 0.03682190383091589, 0.03682190383090726, 0.2994514348233661, -0.004213579037791342, -0.05355060324440208, 0.02318241447626295, 0.02964945857409788, -0.06848925913815337, -0.005389013178829798, 0.3829874114640007, 0.04709386562698295, 0.04709386562697208, 0.3829874114640024, -0.005389013178830043, -0.06848925913815331, 0.02964945857409786, -0.03292028901108043, 0.07604476821551508, 0.005983511330847145, -0.4252373189713157, -0.05228910548955754, -0.05228910548954531, -0.4252373189713177, 0.005983511330847445, 0.07604476821551504, -0.0329202890110804, 0.09489400928644097, -0.219201992388326, -0.01724770337234352, 1.225763056995743, 0.1507253736513569, 0.1507253736513216, 1.225763056995748, -0.01724770337234434, -0.2192019923883258, 0.09489400928644084, -0.4730893437035167, 1.092820584748061, 0.08598756370578201, -6.1109804984601, -0.7514338221809812, -0.7514338221808053, -6.110980498460129, 0.0859875637057862, 1.09282058474806, -0.4730893437035163, 0.4449450582601869, -1.02780822527549, -0.08087212711918088, 5.747435680179285, 0.706730705814681, 0.7067307058145185, 5.747435680179311, -0.08087212711918494, -1.027808225275489, 0.4449450582601864, -0.04925460165060385, 0.1137764849150323, 0.008952396103617317, -0.6362305855168289, -0.0782337925614108, -0.07823379256139273, -0.6362305855168319, 0.008952396103617706, 0.1137764849150322, -0.0492546016506038, 0.002815980818820293, -0.006504821649483379, -0.0005118257962798591, 0.03637453283799384, 0.004472777199565479, 0.004472777199564321, 0.03637453283799404, -0.0005118257962798876, -0.006504821649483374, 0.002815980818820295, -0.004113935537302255, 0.009503053702913197, 0.0007477388759721277, -0.05314044836345182, -0.006534390060027498, -0.006534390060026013, -0.05314044836345205, 0.0007477388759721646, 0.009503053702913192, -0.004113935537302249, 0.002311968049847617, -0.005340568985076732, -0.0004202176663201082, 0.02986410887017532, 0.003672226000393821, 0.003672226000392985, 0.02986410887017545, -0.0004202176663201296, -0.005340568985076729, 0.002311968049847614},
        {-0.02198192343789863, 0.0395907917316884, 0.01136460062590474, 0.02701395152552879, -1.222908076560426, 1.367383721071083, -0.23582853200766, 0.4729675643688158, 0.2042230964705687, -0.03844205209091056, -0.0281140745291355, 0.05063517178361043, 0.01453490773422007, 0.03454985404988617, -1.564054615322062, 1.748833670207433, -0.3016160502828945, 0.6049081824934117, 0.2611938564414026, -0.04916597587968329, 0.0312155264645117, -0.05622107685630401, -0.01613835079531368, -0.03836127994615518, 1.736595959647574, -1.941759230161217, 0.3348893377219305, -0.6716396571992509, -0.2900079008350407, 0.05458980410809167, -0.08997996515788487, 0.1620594335456818, 0.04651942179860217, 0.1105778765862281, -5.005805175833247, 5.597196256598223, -0.9653314985484681, 1.936027349151028, 0.8359590168154002, -0.1573572266098473, 0.4485906221383725, -0.8079392116910494, -0.2319202539091408, -0.5512804807772194, 24.95619168322988, -27.90454237865352, 4.81261197138481, -9.65196765200635, -4.167631925366901, 0.7844965938697586, -0.421903776017874, 0.7598745657688148, 0.2181232197696465, 0.518484571469963, -23.47153726930284, 26.24448933302311, -4.526307646729538, 9.077768008105489, 3.919697736843372, -0.7378265592113, 0.04670397397805059, -0.08411672036029785, -0.02414584025836875, -0.05739529084690629, 2.598256114692661, -2.905216820872467, 0.5010539525025534, -1.004892264371148, -0.4339033483682332, 0.08167604648378231, -0.002670156502692179, 0.004809115557505792, 0.0013804643778097, 0.003281399761630562, -0.148547326263218, 0.1660968633981562, -0.0286462233407178, 0.05745163389044185, 0.02480709345482784, -0.004669577512512197, 0.003900897212498044, -0.007025755027445292, -0.002016754313060389, -0.004793877501311418, 0.2170164371863513, -0.2426549869946819, 0.04184997121544557, -0.08393250293408698, -0.03624129207802923, 0.006821900470529201, -0.002192243811129142, 0.003948365501167157, 0.001133384685761406, 0.002694084901773473, -0.1219598762589321, 0.1363683441271861, -0.02351903559495272, 0.04716876659049449, 0.0203670447944219, -0.003833802397751369},
        {0.03630860018450288, -0.06460783081261466, 0.04422393366747204, -0.7735252391124541, 6.987697007968715, -7.429692565696148, 1.490275197926278, -0.5170009212455957, 0.4656337431327625, 0.364071215683772, 0.04643736907370828, -0.0826310479955554, 0.05656079053370181, -0.9893104342762556, 8.937008402565407, -9.502304523562355, 1.906007365633277, -0.6612252322956808, 0.5955284938871686, 0.4656337431327617, -0.05156018640273594, 0.09174663255667587, -0.0628003903144224, 1.098447897004903, -9.922909680513762, 10.55056740429829, -2.116271809047853, 0.7341694181090087, -0.6612252322956804, -0.5170009212455947, 0.1486242361257834, -0.2644632250616793, 0.1810245596471856, -3.166318646351211, 28.60317183285519, -30.41241955390248, 6.100235530521222, -2.11627180904785, 1.906007365633273, 1.490275197926274, -0.7409587059910442, 1.318468199615025, -0.9024889006343422, 15.78552346583736, -142.5996845532996, 151.6195987015541, -30.41241955390246, 10.55056740429827, -9.502304523562337, -7.429692565696119, 0.696878758723829, -1.240031967956479, 0.848799453668468, -14.84643598858204, 134.1163689182893, -142.5996845532997, 28.60317183285518, -9.922909680513749, 8.937008402565397, 6.987697007968694, -0.07714320009289291, 0.1372692638829018, -0.0939605423086715, 1.643473226577982, -14.84643598858202, 15.78552346583735, -3.166318646351208, 1.098447897004901, -0.9893104342762526, -0.7735252391124511, 0.004410425919287466, -0.007847949250506119, 0.005371905893059423, -0.09396054230867115, 0.8487994536684633, -0.9024889006343368, 0.1810245596471848, -0.06280039031442194, 0.0565607905337014, 0.04422393366747173, -0.006443299543353017, 0.01146526180179283, -0.00784794925050614, 0.1372692638829018, -1.240031967956478, 1.318468199615025, -0.2644632250616792, 0.09174663255667578, -0.08263104799555525, -0.06460783081261443, 0.003621034540954074, -0.006443299543353017, 0.004410425919287493, -0.07714320009289295, 0.696878758723829, -0.7409587059910443, 0.1486242361257835, -0.05156018640273589, 0.04643736907370821, 0.03630860018450276},
        {-0.03844205209091059, -0.04916597587968334, 0.05458980410809173, -0.1573572266098474, 0.7844965938697591, -0.7378265592113006, 0.08167604648378235, -0.0046695775125122, 0.006821900470529211, -0.003833802397751371, 0.2042230964705688, 0.2611938564414029, -0.2900079008350409, 0.8359590168154006, -4.167631925366903, 3.919697736843375, -0.4339033483682336, 0.02480709345482786, -0.03624129207802923, 0.02036704479442191, 0.4729675643688154, 0.6049081824934112, -0.6716396571992499, 1.936027349151026, -9.651967652006341, 9.07776800810548, -1.004892264371147, 0.05745163389044181, -0.08393250293408694, 0.04716876659049443, -0.2358285320076576, -0.3016160502828915, 0.3348893377219272, -0.9653314985484585, 4.812611971384761, -4.526307646729492, 0.5010539525025487, -0.02864622334071753, 0.04184997121544519, -0.02351903559495248, 1.367383721071069, 1.748833670207415, -1.941759230161198, 5.597196256598169, -27.90454237865326, 26.24448933302284, -2.905216820872441, 0.1660968633981547, -0.2426549869946797, 0.1363683441271847, -1.222908076560414, -1.564054615322047, 1.736595959647557, -5.005805175833199, 24.95619168322964, -23.47153726930259, 2.598256114692637, -0.1485473262632165, 0.217016437186349, -0.1219598762589309, 0.02701395152552693, 0.03454985404988387, -0.03836127994615254, 0.1105778765862206, -0.5512804807771818, 0.5184845714699281, -0.05739529084690236, 0.003281399761630343, -0.004793877501311094, 0.002694084901773292, 0.01136460062590503, 0.01453490773422042, -0.01613835079531409, 0.04651942179860329, -0.2319202539091468, 0.2181232197696522, -0.02414584025836937, 0.001380464377809732, -0.002016754313060442, 0.001133384685761436, 0.03959079173168834, 0.05063517178361034, -0.05622107685630394, 0.1620594335456816, -0.8079392116910481, 0.7598745657688142, -0.08411672036029773, 0.004809115557505783, -0.007025755027445284, 0.003948365501167151, -0.02198192343789861, -0.02811407452913547, 0.03121552646451166, -0.08997996515788478, 0.4485906221383719, -0.4219037760178737, 0.04670397397805057, -0.002670156502692174, 0.003900897212498039, -0.00219224381112914},
        {0.004059072250973766, -0.02156378910083279, -0.04994034947978323, 0.02490098728331148, -0.144381192385543, 0.1291260993903018, -0.002852386255736418, -0.001199981076320255, -0.004180366951590879, 0.002321057555376996, -0.02156378910083279, 0.1145574583634581, 0.2653077100427525, -0.1322862972077011, 0.7670239380394658, -0.6859813776405819, 0.01515327933324137, 0.006374889939087098, 0.02220816623467874, -0.01233059983178101, -0.04994034947978324, 0.2653077100427525, 0.6144356030037577, -0.3063665612312971, 1.776378137716584, -1.588688777088381, 0.03509402091293989, 0.01476383533359955, 0.0514326855024818, -0.02855687662381739, 0.02490098728331148, -0.1322862972077011, -0.3063665612312971, 0.1527588397902718, -0.8857280711569646, 0.7921434160454163, -0.01749839113214742, -0.007361463820826728, -0.02564508780945392, 0.01423887555990442, -0.144381192385543, 0.7670239380394656, 1.776378137716583, -0.8857280711569649, 5.135638743476484, -4.59301511412944, 0.1014593737888161, 0.04268332464336941, 0.1486956446599377, -0.0825600129094549, 0.1291260993903018, -0.6859813776405818, -1.588688777088381, 0.7921434160454164, -4.593015114129434, 4.107724256387436, -0.09073933361727983, -0.03817347072110836, -0.1329846933248259, 0.07383684991424355, -0.002852386255736422, 0.01515327933324138, 0.03509402091293989, -0.01749839113214744, 0.1014593737888163, -0.09073933361727995, 0.002004425358519349, 0.000843249224848997, 0.002937622318447773, -0.001631050708235673, -0.001199981076320256, 0.006374889939087104, 0.01476383533359956, -0.007361463820826732, 0.04268332464336945, -0.03817347072110836, 0.0008432492248489987, 0.0003547496803441399, 0.001235839355355215, -0.0006861728422878838, -0.004180366951590879, 0.02220816623467874, 0.0514326855024818, -0.02564508780945392, 0.1486956446599378, -0.1329846933248259, 0.00293762231844778, 0.001235839355355213, 0.004305286225383368, -0.002390416257042631, 0.002321057555376996, -0.01233059983178101, -0.02855687662381739, 0.01423887555990442, -0.08256001290945489, 0.07383684991424358, -0.001631050708235672, -0.0006861728422878836, -0.00239041625704263, 0.001327226479913047},
        {-0.002447816653716587, 0.005654374723238392, 0.0004449091767821289, -0.03161888982230469, -0.003888001808920988, -0.003888001808920082, -0.03161888982230483, 0.0004449091767821508, 0.005654374723238385, -0.002447816653716585, 0.01300400653513568, -0.03003882081668827, -0.002363576468700052, 0.1679751011002768, 0.02065497873587135, 0.02065497873586651, 0.1679751011002775, -0.002363576468700179, -0.03003882081668824, 0.01300400653513567, 0.03011644326353481, -0.06956797817541559, -0.005473891175485211, 0.3890195370407225, 0.04783560309106662, 0.04783560309105555, 0.3890195370407241, -0.00547389117548547, -0.06956797817541552, 0.03011644326353478, -0.0150164982531297, 0.03468760947644118, 0.002729362048340622, -0.1939708201027348, -0.02385153000866334, -0.02385153000865775, -0.1939708201027357, 0.002729362048340757, 0.03468760947644114, -0.01501649825312968, 0.08706883380063168, -0.2011260983442668, -0.01582541858713256, 1.124683851921158, 0.1382962170812537, 0.138296217081221, 1.124683851921163, -0.01582541858713329, -0.2011260983442667, 0.0870688338006316, -0.07786927577877371, 0.1798754265405681, 0.01415332938876472, -1.005851499397787, -0.1236840531448939, -0.1236840531448652, -1.005851499397792, 0.01415332938876544, 0.179875426540568, -0.07786927577877364, 0.001720126705788231, -0.003973435245327163, -0.0003126460290525531, 0.02221918733502457, 0.002732171845274314, 0.002732171845273565, 0.02221918733502473, -0.0003126460290525648, -0.003973435245327145, 0.001720126705788229, 0.0007236465579189474, -0.001671599382021817, -0.0001315282310364643, 0.00934747328824178, 0.001149407624928298, 0.001149407624928024, 0.009347473288241832, -0.0001315282310364748, -0.001671599382021817, 0.0007236465579189459, 0.002520963217714535, -0.005823340843284086, -0.0004582041176116863, 0.03256373720064881, 0.004004184519218607, 0.004004184519217614, 0.03256373720064897, -0.0004582041176117099, -0.005823340843284078, 0.002520963217714533, -0.00139970983195077, 0.003233282967347209, 0.0002544078406042124, -0.01808030470437847, -0.002223236103213126, -0.002223236103212633, -0.01808030470437855, 0.0002544078406042261, 0.003233282967347208, -0.001399709831950768},
        {0.002321057555376999, -0.004180366951590886, -0.001199981076320225, -0.002852386255736611, 0.1291260993903031, -0.1443811923855444, 0.02490098728331172, -0.0499403494797833, -0.02156378910083279, 0.004059072250973763, -0.01233059983178102, 0.02220816623467878, 0.006374889939086937, 0.0151532793332424, -0.6859813776405886, 0.7670239380394731, -0.1322862972077024, 0.2653077100427526, 0.114557458363458, -0.02156378910083279, -0.02855687662381741, 0.05143268550248189, 0.01476383533359917, 0.03509402091294223, -1.588688777088397, 1.776378137716601, -0.3063665612313002, 0.6144356030037583, 0.2653077100427523, -0.0499403494797832, 0.01423887555990444, -0.02564508780945395, -0.007361463820826538, -0.01749839113214859, 0.792143416045424, -0.8857280711569734, 0.1527588397902733, -0.3063665612312975, -0.132286297207701, 0.02490098728331145, -0.082560012909455, 0.1486956446599379, 0.04268332464336837, 0.101459373788823, -4.593015114129484, 5.135638743476537, -0.8857280711569736, 1.776378137716584, 0.767023938039465, -0.1443811923855429, 0.07383684991424366, -0.1329846933248261, -0.03817347072110738, -0.09073933361728594, 4.107724256387476, -4.593015114129482, 0.7921434160454245, -1.588688777088383, -0.6859813776405814, 0.1291260993903017, -0.001631050708235678, 0.002937622318447776, 0.0008432492248489757, 0.002004425358519482, -0.09073933361728079, 0.1014593737888171, -0.01749839113214756, 0.03509402091293994, 0.01515327933324137, -0.00285238625573642, -0.0006861728422878832, 0.001235839355355216, 0.0003547496803441341, 0.0008432492248490559, -0.03817347072110878, 0.04268332464336991, -0.007361463820826802, 0.01476383533359957, 0.006374889939087101, -0.001199981076320254, -0.002390416257042633, 0.004305286225383374, 0.001235839355355182, 0.002937622318447972, -0.1329846933248271, 0.1486956446599391, -0.02564508780945419, 0.05143268550248186, 0.02220816623467873, -0.004180366951590877, 0.001327226479913048, -0.002390416257042635, -0.0006861728422878665, -0.001631050708235782, 0.07383684991424433, -0.08256001290945571, 0.01423887555990457, -0.02855687662381742, -0.01233059983178101, 0.002321057555376995},
        {-0.003833802397751383, 0.006821900470529232, -0.00466957751251224, 0.08167604648378267, -0.7378265592113027, 0.7844965938697619, -0.1573572266098478, 0.05458980410809181, -0.04916597587968339, -0.03844205209091058, 0.02036704479442197, -0.03624129207802936, 0.02480709345482805, -0.4339033483682351, 3.919697736843385, -4.167631925366917, 0.8359590168154029, -0.2900079008350414, 0.261193856441403, 0.2042230964705688, 0.04716876659049458, -0.08393250293408719, 0.05745163389044231, -1.004892264371151, 9.077768008105505, -9.651967652006377, 1.93602734915103, -0.6716396571992511, 0.6049081824934118, 0.4729675643688152, -0.02351903559495256, 0.0418499712154453, -0.02864622334071778, 0.5010539525025504, -4.526307646729505, 4.812611971384777, -0.9653314985484608, 0.3348893377219276, -0.3016160502828917, -0.2358285320076575, 0.1363683441271852, -0.2426549869946804, 0.166096863398156, -2.90521682087245, 26.24448933302293, -27.90454237865336, 5.597196256598183, -1.941759230161201, 1.748833670207418, 1.367383721071068, -0.1219598762589313, 0.2170164371863498, -0.1485473262632178, 2.598256114692647, -23.47153726930266, 24.95619168322973, -5.005805175833212, 1.73659595964756, -1.564054615322049, -1.222908076560414, 0.002694084901773293, -0.004793877501311118, 0.003281399761630367, -0.05739529084690262, 0.5184845714699289, -0.5512804807771841, 0.1105778765862209, -0.03836127994615264, 0.03454985404988391, 0.02701395152552688, 0.001133384685761441, -0.002016754313060447, 0.00138046437780974, -0.02414584025836942, 0.2181232197696525, -0.2319202539091475, 0.04651942179860345, -0.01613835079531412, 0.01453490773422042, 0.01136460062590501, 0.003948365501167165, -0.007025755027445305, 0.004809115557505817, -0.08411672036029803, 0.7598745657688155, -0.8079392116910508, 0.1620594335456819, -0.05622107685630402, 0.05063517178361041, 0.03959079173168832, -0.002192243811129146, 0.003900897212498052, -0.0026701565026922, 0.04670397397805074, -0.4219037760178749, 0.4485906221383735, -0.08997996515788503, 0.03121552646451171, -0.02811407452913552, -0.0219819234378986},
        {0.02318241447626297, 0.02964945857409789, -0.03292028901108042, 0.09489400928644098, -0.4730893437035169, 0.4449450582601869, -0.04925460165060386, 0.002815980818820296, -0.004113935537302255, 0.002311968049847618, -0.05355060324440213, -0.06848925913815339, 0.07604476821551509, -0.219201992388326, 1.092820584748061, -1.02780822527549, 0.1137764849150324, -0.006504821649483383, 0.009503053702913199, -0.005340568985076732, -0.004213579037791124, -0.00538901317882978, 0.005983511330847166, -0.01724770337234351, 0.08598756370578226, -0.0808721271191811, 0.008952396103617314, -0.0005118257962798601, 0.0007477388759721251, -0.0004202176663201088, 0.2994514348233648, 0.3829874114640006, -0.4252373189713158, 1.225763056995744, -6.1109804984601, 5.747435680179287, -0.636230585516829, 0.03637453283799389, -0.05314044836345183, 0.02986410887017531, 0.03682190383091568, 0.04709386562698289, -0.05228910548955704, 0.150725373651356, -0.7514338221809812, 0.7067307058146766, -0.07823379256141041, 0.004472777199565396, -0.006534390060027456, 0.003672226000393821, 0.03682190383090746, 0.04709386562697222, -0.05228910548954555, 0.150725373651323, -0.7514338221808113, 0.7067307058145229, -0.07823379256139287, 0.004472777199564425, -0.006534390060026027, 0.003672226000392978, 0.2994514348233662, 0.3829874114640024, -0.4252373189713177, 1.225763056995748, -6.110980498460129, 5.747435680179311, -0.6362305855168318, 0.03637453283799402, -0.05314044836345205, 0.02986410887017545, -0.004213579037791331, -0.005389013178830044, 0.005983511330847424, -0.01724770337234431, 0.08598756370578606, -0.08087212711918476, 0.008952396103617661, -0.0005118257962798837, 0.0007477388759721546, -0.0004202176663201328, -0.05355060324440208, -0.06848925913815333, 0.07604476821551504, -0.2192019923883258, 1.092820584748061, -1.027808225275489, 0.1137764849150323, -0.006504821649483378, 0.009503053702913194, -0.005340568985076727, 0.02318241447626294, 0.02964945857409785, -0.03292028901108039, 0.09489400928644087, -0.4730893437035161, 0.4449450582601863, -0.0492546016506038, 0.002815980818820292, -0.004113935537302247, 0.002311968049847615},
        {-0.002447816653716587, 0.01300400653513568, 0.03011644326353481, -0.0150164982531297, 0.08706883380063168, -0.07786927577877369, 0.001720126705788228, 0.0007236465579189473, 0.002520963217714535, -0.00139970983195077, 0.005654374723238391, -0.03003882081668827, -0.06956797817541559, 0.03468760947644118, -0.2011260983442668, 0.1798754265405682, -0.003973435245327161, -0.001671599382021817, -0.005823340843284087, 0.00323328296734721, 0.0004449091767821277, -0.002363576468700042, -0.005473891175485174, 0.002729362048340625, -0.01582541858713248, 0.0141533293887647, -0.0003126460290525542, -0.0001315282310364657, -0.0004582041176116849, 0.0002544078406042123, -0.03161888982230467, 0.1679751011002768, 0.3890195370407225, -0.1939708201027348, 1.124683851921158, -1.005851499397788, 0.02221918733502459, 0.009347473288241787, 0.03256373720064882, -0.01808030470437846, -0.003888001808920989, 0.0206549787358714, 0.04783560309106695, -0.02385153000866328, 0.138296217081253, -0.1236840531448935, 0.002732171845274203, 0.001149407624928316, 0.004004184519218579, -0.002223236103213153, -0.003888001808920074, 0.02065497873586655, 0.04783560309105528, -0.02385153000865786, 0.1382962170812221, -0.1236840531448658, 0.00273217184527369, 0.001149407624927982, 0.004004184519217649, -0.002223236103212609, -0.03161888982230483, 0.1679751011002776, 0.3890195370407241, -0.1939708201027357, 1.124683851921163, -1.005851499397792, 0.02221918733502469, 0.009347473288241837, 0.03256373720064895, -0.01808030470437855, 0.0004449091767821522, -0.002363576468700179, -0.005473891175485476, 0.002729362048340763, -0.01582541858713327, 0.01415332938876539, -0.000312646029052557, -0.0001315282310364759, -0.0004582041176117069, 0.0002544078406042275, 0.005654374723238388, -0.03003882081668825, -0.0695679781754155, 0.03468760947644115, -0.2011260983442666, 0.1798754265405681, -0.003973435245327164, -0.001671599382021813, -0.005823340843284082, 0.003233282967347208, -0.002447816653716585, 0.01300400653513567, 0.03011644326353477, -0.01501649825312968, 0.0870688338006316, -0.07786927577877364, 0.001720126705788225, 0.0007236465579189479, 0.002520963217714531, -0.001399709831950769},
        {0.00147615169175046, -0.003409861110645825, -0.0002683017263507356, 0.0190677180137653, 0.002344652913057627, 0.002344652913057076, 0.01906771801376538, -0.0002683017263507509, -0.003409861110645824, 0.001476151691750458, -0.003409861110645825, 0.007876665290480407, 0.0006197680277137592, -0.04404579182969925, -0.00541606992755402, -0.005416069927552732, -0.04404579182969945, 0.0006197680277137934, 0.0078766652904804, -0.003409861110645821, -0.0002683017263507354, 0.0006197680277137599, 4.87658665197341e-05, -0.003465701857913854, -0.0004261583872323359, -0.0004261583872322127, -0.003465701857913857, 4.876586651973789e-05, 0.0006197680277137609, -0.0002683017263507361, 0.0190677180137653, -0.04404579182969924, -0.003465701857913862, 0.2463011574517314, 0.03028630515155337, 0.0302863051515463, 0.2463011574517324, -0.003465701857914041, -0.04404579182969923, 0.01906771801376528, 0.002344652913057636, -0.005416069927554032, -0.0004261583872323116, 0.03028630515155329, 0.003724141166136452, 0.003724141166135286, 0.03028630515155342, -0.0004261583872323393, -0.005416069927553994, 0.002344652913057629, 0.002344652913057091, -0.005416069927552752, -0.0004261583872322509, 0.03028630515154637, 0.003724141166135564, 0.003724141166134509, 0.03028630515154633, -0.0004261583872322283, -0.005416069927552766, 0.002344652913057074, 0.01906771801376538, -0.04404579182969945, -0.003465701857913868, 0.2463011574517324, 0.03028630515155346, 0.03028630515154644, 0.2463011574517335, -0.003465701857914058, -0.04404579182969942, 0.01906771801376536, -0.0002683017263507514, 0.0006197680277137933, 4.876586651973594e-05, -0.003465701857914037, -0.0004261583872323185, -0.0004261583872322491, -0.003465701857914059, 4.876586651974386e-05, 0.0006197680277137904, -0.0002683017263507505, -0.003409861110645823, 0.0078766652904804, 0.0006197680277137621, -0.04404579182969923, -0.005416069927553988, -0.005416069927552721, -0.04404579182969942, 0.0006197680277137904, 0.007876665290480396, -0.003409861110645819, 0.001476151691750459, -0.003409861110645821, -0.0002683017263507349, 0.01906771801376528, 0.002344652913057618, 0.002344652913057077, 0.01906771801376536, -0.0002683017263507509, -0.00340986111064582, 0.001476151691750457},
        {-0.001399709831950772, 0.002520963217714538, 0.000723646557918929, 0.001720126705788342, -0.07786927577877449, 0.08706883380063253, -0.01501649825312985, 0.03011644326353485, 0.01300400653513568, -0.002447816653716586, 0.003233282967347214, -0.005823340843284092, -0.001671599382021772, -0.003973435245327424, 0.17987542654057, -0.2011260983442687, 0.03468760947644151, -0.06956797817541563, -0.03003882081668826, 0.005654374723238387, 0.0002544078406042138, -0.0004582041176116853, -0.000131528231036462, -0.0003126460290525672, 0.01415332938876489, -0.01582541858713267, 0.002729362048340641, -0.005473891175485221, -0.002363576468700056, 0.0004449091767821273, -0.01808030470437849, 0.03256373720064887, 0.009347473288241559, 0.02221918733502606, -1.005851499397797, 1.124683851921169, -0.1939708201027368, 0.3890195370407228, 0.1679751011002767, -0.03161888982230465, -0.002223236103213143, 0.004004184519218601, 0.001149407624928243, 0.002732171845274523, -0.1236840531448961, 0.1382962170812553, -0.02385153000866352, 0.0478356030910667, 0.0206549787358713, -0.003888001808921019, -0.00222323610321261, 0.004004184519217619, 0.001149407624928031, 0.002732171845273787, -0.1236840531448661, 0.1382962170812216, -0.02385153000865792, 0.04783560309105531, 0.02065497873586651, -0.003888001808920044, -0.01808030470437857, 0.03256373720064901, 0.009347473288241582, 0.02221918733502619, -1.005851499397802, 1.124683851921174, -0.1939708201027377, 0.3890195370407246, 0.1679751011002776, -0.0316188898223048, 0.000254407840604227, -0.0004582041176117106, -0.0001315282310364678, -0.0003126460290525874, 0.01415332938876562, -0.0158254185871335, 0.002729362048340801, -0.005473891175485473, -0.002363576468700171, 0.0004449091767821526, 0.003233282967347212, -0.005823340843284087, -0.001671599382021776, -0.003973435245327415, 0.1798754265405698, -0.2011260983442685, 0.03468760947644147, -0.06956797817541556, -0.03003882081668824, 0.005654374723238382, -0.00139970983195077, 0.002520963217714537, 0.0007236465579189274, 0.001720126705788345, -0.07786927577877444, 0.0870688338006325, -0.01501649825312984, 0.03011644326353482, 0.01300400653513567, -0.002447816653716585},
        {0.002311968049847624, -0.004113935537302263, 0.002815980818820318, -0.04925460165060401, 0.444945058260188, -0.4730893437035182, 0.09489400928644122, -0.03292028901108048, 0.02964945857409792, 0.02318241447626295, -0.005340568985076747, 0.009503053702913223, -0.006504821649483435, 0.1137764849150328, -1.027808225275492, 1.092820584748065, -0.2192019923883265, 0.07604476821551522, -0.06848925913815346, -0.05355060324440211, -0.0004202176663201115, 0.0007477388759721225, -0.0005118257962798589, 0.008952396103617307, -0.08087212711918085, 0.08598756370578201, -0.01724770337234349, 0.005983511330847124, -0.005389013178829774, -0.004213579037791124, 0.02986410887017541, -0.05314044836345196, 0.03637453283799416, -0.6362305855168311, 5.747435680179304, -6.110980498460117, 1.225763056995746, -0.4252373189713166, 0.3829874114640009, 0.2994514348233646, 0.003672226000393842, -0.006534390060027581, 0.004472777199565542, -0.0782337925614115, 0.7067307058146888, -0.7514338221809869, 0.1507253736513583, -0.05228910548955779, 0.04709386562698348, 0.0368219038309161, 0.003672226000392965, -0.006534390060025951, 0.004472777199564348, -0.07823379256139251, 0.7067307058145129, -0.7514338221808039, 0.1507253736513208, -0.05228910548954491, 0.04709386562697182, 0.03682190383090705, 0.02986410887017554, -0.05314044836345222, 0.03637453283799436, -0.636230585516834, 5.747435680179327, -6.110980498460148, 1.225763056995752, -0.4252373189713183, 0.3829874114640028, 0.299451434823366, -0.0004202176663201333, 0.0007477388759721654, -0.0005118257962798915, 0.008952396103617748, -0.08087212711918507, 0.08598756370578658, -0.01724770337234441, 0.005983511330847475, -0.00538901317883007, -0.004213579037791359, -0.005340568985076744, 0.009503053702913211, -0.006504821649483428, 0.1137764849150326, -1.027808225275491, 1.092820584748063, -0.2192019923883262, 0.07604476821551512, -0.06848925913815337, -0.05355060324440204, 0.002311968049847621, -0.004113935537302263, 0.002815980818820321, -0.049254601650604, 0.4449450582601877, -0.4730893437035179, 0.09489400928644118, -0.03292028901108045, 0.0296494585740979, 0.02318241447626294},
        {-0.02198192343789863, -0.02811407452913551, 0.03121552646451169, -0.08997996515788485, 0.4485906221383725, -0.421903776017874, 0.04670397397805059, -0.002670156502692178, 0.003900897212498043, -0.002192243811129142, 0.03959079173168839, 0.05063517178361043, -0.05622107685630401, 0.1620594335456818, -0.8079392116910495, 0.7598745657688148, -0.08411672036029784, 0.004809115557505789, -0.007025755027445295, 0.003948365501167157, 0.01136460062590474, 0.01453490773422008, -0.01613835079531364, 0.04651942179860211, -0.231920253909141, 0.2181232197696466, -0.0241458402583688, 0.001380464377809692, -0.002016754313060393, 0.001133384685761408, 0.0270139515255288, 0.03454985404988616, -0.03836127994615519, 0.1105778765862281, -0.5512804807772196, 0.5184845714699627, -0.05739529084690638, 0.003281399761630593, -0.004793877501311421, 0.002694084901773468, -1.222908076560426, -1.564054615322062, 1.736595959647574, -5.005805175833247, 24.95619168322988, -23.47153726930283, 2.598256114692662, -0.1485473262632181, 0.2170164371863512, -0.1219598762589321, 1.367383721071082, 1.748833670207433, -1.941759230161217, 5.597196256598224, -27.90454237865352, 26.2444893330231, -2.905216820872468, 0.1660968633981564, -0.242654986994682, 0.136368344127186, -0.23582853200766, -0.3016160502828945, 0.3348893377219305, -0.9653314985484683, 4.81261197138481, -4.526307646729538, 0.5010539525025535, -0.02864622334071784, 0.0418499712154456, -0.02351903559495272, 0.4729675643688159, 0.6049081824934114, -0.6716396571992509, 1.936027349151027, -9.651967652006348, 9.077768008105485, -1.004892264371148, 0.0574516338904419, -0.08393250293408699, 0.04716876659049449, 0.2042230964705687, 0.2611938564414027, -0.2900079008350407, 0.8359590168154001, -4.167631925366902, 3.919697736843373, -0.4339033483682332, 0.02480709345482781, -0.03624129207802922, 0.02036704479442191, -0.03844205209091056, -0.04916597587968329, 0.05458980410809169, -0.1573572266098472, 0.7844965938697587, -0.7378265592113, 0.08167604648378235, -0.004669577512512206, 0.006821900470529201, -0.003833802397751367},
        {0.002321057555376998, -0.01233059983178102, -0.02855687662381741, 0.01423887555990444, -0.08256001290945497, 0.07383684991424366, -0.001631050708235671, -0.0006861728422878842, -0.002390416257042633, 0.001327226479913048, -0.004180366951590887, 0.02220816623467878, 0.05143268550248189, -0.02564508780945396, 0.148695644659938, -0.132984693324826, 0.002937622318447778, 0.001235839355355215, 0.004305286225383375, -0.002390416257042635, -0.001199981076320223, 0.006374889939086958, 0.01476383533359917, -0.007361463820826536, 0.04268332464336842, -0.03817347072110737, 0.0008432492248489705, 0.0003547496803441341, 0.001235839355355181, -0.0006861728422878674, -0.002852386255736616, 0.01515327933324241, 0.03509402091294225, -0.0174983911321486, 0.1014593737888232, -0.09073933361728607, 0.002004425358519496, 0.0008432492248490429, 0.002937622318447974, -0.001631050708235781, 0.1291260993903031, -0.6859813776405884, -1.588688777088397, 0.792143416045424, -4.593015114129484, 4.107724256387477, -0.09073933361728079, -0.03817347072110869, -0.1329846933248271, 0.07383684991424433, -0.1443811923855444, 0.7670239380394731, 1.776378137716601, -0.8857280711569735, 5.135638743476537, -4.593015114129482, 0.1014593737888172, 0.04268332464336989, 0.1486956446599391, -0.08256001290945572, 0.02490098728331172, -0.1322862972077024, -0.3063665612313001, 0.1527588397902733, -0.8857280711569737, 0.7921434160454242, -0.01749839113214759, -0.007361463820826791, -0.02564508780945418, 0.01423887555990457, -0.0499403494797833, 0.2653077100427527, 0.6144356030037584, -0.3063665612312974, 1.776378137716585, -1.588688777088383, 0.03509402091293988, 0.01476383533359957, 0.05143268550248186, -0.02855687662381742, -0.02156378910083279, 0.114557458363458, 0.2653077100427523, -0.132286297207701, 0.7670239380394649, -0.6859813776405816, 0.01515327933324136, 0.006374889939087109, 0.02220816623467873, -0.01233059983178101, 0.004059072250973763, -0.02156378910083278, -0.0499403494797832, 0.02490098728331145, -0.1443811923855429, 0.1291260993903017, -0.002852386255736417, -0.001199981076320253, -0.004180366951590876, 0.002321057555376994},
        {-0.001399709831950772, 0.003233282967347214, 0.0002544078406042144, -0.01808030470437848, -0.002223236103213126, -0.002223236103212627, -0.01808030470437857, 0.000254407840604227, 0.003233282967347211, -0.00139970983195077, 0.002520963217714538, -0.005823340843284092, -0.000458204117611684, 0.03256373720064886, 0.004004184519218621, 0.004004184519217627, 0.032563737200649, -0.0004582041176117117, -0.005823340843284088, 0.002520963217714536, 0.0007236465579189285, -0.001671599382021772, -0.0001315282310364612, 0.009347473288241559, 0.001149407624928244, 0.001149407624928019, 0.009347473288241577, -0.0001315282310364683, -0.001671599382021776, 0.0007236465579189269, 0.001720126705788341, -0.003973435245327419, -0.0003126460290525694, 0.02221918733502608, 0.002732171845274425, 0.002732171845273842, 0.02221918733502619, -0.0003126460290525839, -0.003973435245327418, 0.001720126705788341, -0.07786927577877449, 0.17987542654057, 0.01415332938876484, -1.005851499397797, -0.1236840531448959, -0.1236840531448666, -1.005851499397802, 0.01415332938876556, 0.1798754265405698, -0.07786927577877442, 0.08706883380063253, -0.2011260983442687, -0.0158254185871327, 1.124683851921169, 0.1382962170812551, 0.1382962170812223, 1.124683851921174, -0.01582541858713345, -0.2011260983442685, 0.08706883380063248, -0.01501649825312985, 0.03468760947644151, 0.002729362048340651, -0.1939708201027368, -0.02385153000866362, -0.023851530008658, -0.1939708201027376, 0.002729362048340789, 0.03468760947644148, -0.01501649825312984, 0.03011644326353485, -0.06956797817541563, -0.005473891175485217, 0.3890195370407228, 0.04783560309106659, 0.04783560309105545, 0.3890195370407245, -0.0054738911754855, -0.06956797817541556, 0.03011644326353482, 0.01300400653513569, -0.03003882081668825, -0.002363576468700059, 0.1679751011002768, 0.02065497873587141, 0.02065497873586645, 0.1679751011002776, -0.002363576468700177, -0.03003882081668822, 0.01300400653513567, -0.002447816653716586, 0.005654374723238386, 0.0004449091767821289, -0.03161888982230466, -0.003888001808920981, -0.003888001808920086, -0.0316188898223048, 0.0004449091767821508, 0.005654374723238383, -0.002447816653716584},
        {0.001327226479913049, -0.002390416257042635, -0.000686172842287868, -0.001631050708235781, 0.07383684991424436, -0.08256001290945578, 0.01423887555990457, -0.02855687662381744, -0.01233059983178102, 0.002321057555376995, -0.002390416257042636, 0.004305286225383383, 0.001235839355355182, 0.002937622318447979, -0.1329846933248273, 0.1486956446599394, -0.02564508780945423, 0.05143268550248194, 0.02220816623467876, -0.004180366951590884, -0.0006861728422878671, 0.001235839355355185, 0.0003547496803441227, 0.0008432492248490356, -0.03817347072110777, 0.04268332464336891, -0.007361463820826633, 0.01476383533359923, 0.006374889939086961, -0.001199981076320225, -0.001631050708235779, 0.002937622318447975, 0.0008432492248490299, 0.002004425358519623, -0.09073933361728681, 0.1014593737888241, -0.01749839113214877, 0.03509402091294216, 0.01515327933324235, -0.002852386255736609, 0.07383684991424436, -0.1329846933248273, -0.03817347072110774, -0.09073933361728689, 4.107724256387515, -4.593015114129527, 0.7921434160454319, -1.588688777088398, -0.6859813776405878, 0.1291260993903029, -0.08256001290945576, 0.1486956446599394, 0.04268332464336884, 0.1014593737888239, -4.593015114129527, 5.135638743476585, -0.8857280711569819, 1.776378137716601, 0.7670239380394723, -0.1443811923855443, 0.01423887555990457, -0.02564508780945422, -0.007361463820826622, -0.01749839113214878, 0.7921434160454319, -0.8857280711569822, 0.1527588397902748, -0.3063665612313003, -0.1322862972077023, 0.02490098728331169, -0.02855687662381743, 0.05143268550248194, 0.01476383533359921, 0.03509402091294223, -1.588688777088398, 1.776378137716602, -0.3063665612313004, 0.6144356030037588, 0.2653077100427526, -0.04994034947978326, -0.01233059983178102, 0.02220816623467876, 0.00637488993908695, 0.01515327933324238, -0.685981377640588, 0.7670239380394727, -0.1322862972077024, 0.2653077100427526, 0.1145574583634579, -0.02156378910083278, 0.002321057555376995, -0.004180366951590883, -0.001199981076320224, -0.002852386255736611, 0.1291260993903029, -0.1443811923855443, 0.02490098728331169, -0.04994034947978324, -0.02156378910083277, 0.004059072250973759},
        {-0.002192243811129149, 0.003900897212498052, -0.002670156502692199, 0.04670397397805075, -0.4219037760178749, 0.4485906221383736, -0.08997996515788503, 0.03121552646451173, -0.02811407452913552, -0.02198192343789861, 0.003948365501167169, -0.007025755027445318, 0.004809115557505832, -0.08411672036029816, 0.7598745657688173, -0.8079392116910523, 0.1620594335456822, -0.05622107685630412, 0.05063517178361047, 0.03959079173168838, 0.001133384685761408, -0.002016754313060394, 0.001380464377809711, -0.02414584025836885, 0.2181232197696466, -0.2319202539091413, 0.04651942179860221, -0.0161383507953137, 0.01453490773422007, 0.01136460062590473, 0.00269408490177348, -0.004793877501311446, 0.003281399761630607, -0.05739529084690662, 0.5184845714699649, -0.5512804807772218, 0.1105778765862286, -0.03836127994615529, 0.03454985404988629, 0.02701395152552882, -0.1219598762589325, 0.2170164371863519, -0.1485473262632193, 2.598256114692671, -23.4715372693029, 24.95619168322997, -5.005805175833261, 1.736595959647576, -1.564054615322064, -1.222908076560426, 0.1363683441271865, -0.2426549869946827, 0.1660968633981575, -2.905216820872477, 26.24448933302317, -27.90454237865362, 5.597196256598238, -1.941759230161219, 1.748833670207434, 1.367383721071082, -0.02351903559495279, 0.04184997121544572, -0.02864622334071804, 0.5010539525025552, -4.52630764672955, 4.812611971384825, -0.9653314985484711, 0.334889337721931, -0.3016160502828949, -0.2358285320076599, 0.04716876659049465, -0.08393250293408726, 0.05745163389044232, -1.004892264371151, 9.07776800810551, -9.65196765200638, 1.936027349151032, -0.6716396571992517, 0.6049081824934119, 0.4729675643688157, 0.02036704479442198, -0.03624129207802933, 0.02480709345482799, -0.4339033483682346, 3.91969773684338, -4.167631925366915, 0.8359590168154016, -0.2900079008350412, 0.2611938564414029, 0.2042230964705686, -0.003833802397751379, 0.006821900470529224, -0.00466957751251224, 0.08167604648378263, -0.7378265592113019, 0.7844965938697611, -0.1573572266098477, 0.05458980410809176, -0.04916597587968335, -0.03844205209091055},
        {0.03630860018450288, 0.04643736907370827, -0.05156018640273594, 0.1486242361257834, -0.7409587059910442, 0.696878758723829, -0.0771432000928929, 0.004410425919287468, -0.006443299543353017, 0.003621034540954074, -0.06460783081261467, -0.08263104799555542, 0.09174663255667588, -0.2644632250616793, 1.318468199615025, -1.240031967956479, 0.1372692638829019, -0.007847949250506111, 0.01146526180179283, -0.006443299543353017, 0.04422393366747204, 0.05656079053370181, -0.06280039031442231, 0.1810245596471857, -0.9024889006343425, 0.8487994536684684, -0.0939605423086716, 0.005371905893059416, -0.00784794925050615, 0.004410425919287492, -0.7735252391124541, -0.9893104342762558, 1.098447897004903, -3.166318646351212, 15.78552346583736, -14.84643598858204, 1.643473226577982, -0.09396054230867112, 0.1372692638829019, -0.07714320009289294, 6.987697007968715, 8.937008402565413, -9.922909680513762, 28.60317183285519, -142.5996845532997, 134.1163689182893, -14.84643598858203, 0.8487994536684638, -1.240031967956478, 0.6968787587238289, -7.429692565696147, -9.502304523562353, 10.55056740429829, -30.41241955390248, 151.619598701554, -142.5996845532998, 15.78552346583736, -0.9024889006343375, 1.318468199615025, -0.7409587059910441, 1.490275197926278, 1.906007365633277, -2.116271809047853, 6.100235530521223, -30.41241955390247, 28.60317183285518, -3.16631864635121, 0.1810245596471848, -0.2644632250616792, 0.1486242361257835, -0.5170009212455956, -0.6612252322956809, 0.7341694181090086, -2.11627180904785, 10.55056740429827, -9.922909680513747, 1.098447897004901, -0.06280039031442201, 0.09174663255667576, -0.05156018640273588, 0.4656337431327624, 0.5955284938871686, -0.6612252322956805, 1.906007365633273, -9.502304523562334, 8.937008402565397, -0.9893104342762534, 0.05656079053370145, -0.08263104799555523, 0.0464373690737082, 0.3640712156837718, 0.4656337431327617, -0.5170009212455946, 1.490275197926274, -7.429692565696119, 6.987697007968695, -0.7735252391124511, 0.0442239336674717, -0.06460783081261444, 0.03630860018450276},
        {-0.003833802397751383, 0.02036704479442198, 0.04716876659049459, -0.02351903559495256, 0.1363683441271852, -0.1219598762589313, 0.002694084901773301, 0.001133384685761439, 0.003948365501167164, -0.002192243811129146, 0.006821900470529231, -0.03624129207802936, -0.08393250293408716, 0.04184997121544529, -0.2426549869946804, 0.2170164371863497, -0.004793877501311121, -0.002016754313060446, -0.007025755027445304, 0.003900897212498051, -0.00466957751251224, 0.02480709345482805, 0.05745163389044227, -0.02864622334071779, 0.166096863398156, -0.1485473262632178, 0.003281399761630371, 0.001380464377809741, 0.00480911555750582, -0.002670156502692198, 0.08167604648378268, -0.4339033483682351, -1.00489226437115, 0.5010539525025502, -2.905216820872451, 2.598256114692647, -0.05739529084690259, -0.02414584025836943, -0.08411672036029799, 0.04670397397805073, -0.7378265592113026, 3.919697736843385, 9.077768008105501, -4.526307646729506, 26.24448933302292, -23.47153726930267, 0.5184845714699295, 0.2181232197696522, 0.7598745657688161, -0.4219037760178748, 0.7844965938697617, -4.167631925366918, -9.651967652006375, 4.812611971384777, -27.90454237865335, 24.95619168322973, -0.5512804807771847, -0.2319202539091476, -0.8079392116910511, 0.4485906221383734, -0.1573572266098478, 0.8359590168154027, 1.936027349151031, -0.965331498548461, 5.597196256598185, -5.005805175833212, 0.1105778765862211, 0.04651942179860342, 0.162059433545682, -0.08997996515788502, 0.05458980410809182, -0.2900079008350414, -0.671639657199251, 0.3348893377219276, -1.941759230161201, 1.73659595964756, -0.0383612799461527, -0.01613835079531408, -0.05622107685630403, 0.0312155264645117, -0.04916597587968338, 0.2611938564414031, 0.6049081824934116, -0.3016160502828918, 1.748833670207417, -1.564054615322049, 0.03454985404988389, 0.01453490773422043, 0.05063517178361041, -0.02811407452913551, -0.03844205209091058, 0.2042230964705689, 0.4729675643688152, -0.2358285320076575, 1.367383721071068, -1.222908076560414, 0.02701395152552696, 0.01136460062590501, 0.03959079173168832, -0.0219819234378986},
        {0.002311968049847625, -0.005340568985076748, -0.0004202176663201111, 0.02986410887017541, 0.003672226000393814, 0.003672226000393019, 0.02986410887017554, -0.000420217666320132, -0.005340568985076747, 0.002311968049847621, -0.004113935537302264, 0.009503053702913218, 0.0007477388759721236, -0.05314044836345196, -0.006534390060027546, -0.00653439006002602, -0.05314044836345219, 0.0007477388759721659, 0.009503053702913211, -0.004113935537302261, 0.002815980818820315, -0.006504821649483437, -0.0005118257962798568, 0.03637453283799415, 0.004472777199565486, 0.004472777199564404, 0.03637453283799434, -0.0005118257962798949, -0.006504821649483433, 0.002815980818820316, -0.04925460165060401, 0.1137764849150328, 0.008952396103617272, -0.6362305855168311, -0.07823379256141116, -0.07823379256139273, -0.6362305855168338, 0.008952396103617685, 0.1137764849150327, -0.04925460165060396, 0.4449450582601879, -1.027808225275492, -0.0808721271191809, 5.747435680179303, 0.7067307058146834, 0.70673070581452, 5.747435680179327, -0.08087212711918501, -1.027808225275491, 0.4449450582601874, -0.4730893437035182, 1.092820584748065, 0.08598756370578203, -6.110980498460119, -0.7514338221809833, -0.7514338221808075, -6.110980498460147, 0.08598756370578664, 1.092820584748064, -0.4730893437035177, 0.09489400928644125, -0.2192019923883265, -0.01724770337234344, 1.225763056995746, 0.150725373651357, 0.1507253736513225, 1.225763056995751, -0.01724770337234435, -0.2192019923883263, 0.09489400928644112, -0.03292028901108047, 0.0760447682155152, 0.005983511330847123, -0.4252373189713164, -0.05228910548955745, -0.05228910548954535, -0.4252373189713183, 0.005983511330847447, 0.07604476821551515, -0.03292028901108043, 0.02964945857409792, -0.06848925913815344, -0.005389013178829782, 0.3829874114640011, 0.04709386562698314, 0.04709386562697249, 0.3829874114640027, -0.00538901317883006, -0.0684892591381534, 0.02964945857409788, 0.02318241447626295, -0.0535506032444021, -0.004213579037791119, 0.2994514348233647, 0.03682190383091577, 0.03682190383090739, 0.299451434823366, -0.004213579037791338, -0.05355060324440206, 0.02318241447626293},
        {-0.002192243811129148, 0.00394836550116717, 0.001133384685761411, 0.002694084901773482, -0.1219598762589325, 0.1363683441271865, -0.02351903559495279, 0.04716876659049464, 0.02036704479442197, -0.003833802397751379, 0.003900897212498053, -0.007025755027445315, -0.002016754313060392, -0.004793877501311442, 0.2170164371863519, -0.2426549869946827, 0.04184997121544572, -0.08393250293408726, -0.03624129207802933, 0.006821900470529227, -0.002670156502692197, 0.004809115557505832, 0.001380464377809707, 0.003281399761630596, -0.1485473262632193, 0.1660968633981576, -0.02864622334071803, 0.05745163389044234, 0.02480709345482801, -0.004669577512512235, 0.04670397397805075, -0.08411672036029816, -0.02414584025836883, -0.05739529084690653, 2.598256114692671, -2.905216820872477, 0.5010539525025552, -1.004892264371151, -0.4339033483682347, 0.08167604648378263, -0.4219037760178749, 0.7598745657688172, 0.218123219769647, 0.5184845714699646, -23.4715372693029, 26.24448933302317, -4.52630764672955, 9.077768008105506, 3.919697736843381, -0.7378265592113019, 0.4485906221383736, -0.8079392116910522, -0.2319202539091415, -0.5512804807772222, 24.95619168322997, -27.90454237865362, 4.812611971384825, -9.65196765200638, -4.167631925366914, 0.7844965938697611, -0.08997996515788505, 0.1620594335456822, 0.04651942179860227, 0.1105778765862285, -5.005805175833261, 5.597196256598238, -0.965331498548471, 1.936027349151032, 0.8359590168154019, -0.1573572266098477, 0.03121552646451172, -0.05622107685630412, -0.01613835079531369, -0.03836127994615531, 1.736595959647577, -1.941759230161219, 0.334889337721931, -0.6716396571992516, -0.2900079008350411, 0.05458980410809176, -0.02811407452913552, 0.05063517178361047, 0.01453490773422007, 0.03454985404988628, -1.564054615322064, 1.748833670207434, -0.3016160502828949, 0.6049081824934122, 0.2611938564414029, -0.04916597587968335, -0.02198192343789861, 0.03959079173168838, 0.01136460062590474, 0.0270139515255288, -1.222908076560426, 1.367383721071082, -0.2358285320076599, 0.4729675643688157, 0.2042230964705686, -0.03844205209091055},
        {0.003621034540954083, -0.006443299543353039, 0.00441042591928751, -0.07714320009289317, 0.6968787587238311, -0.7409587059910465, 0.1486242361257839, -0.05156018640273603, 0.04643736907370834, 0.03630860018450287, -0.006443299543353038, 0.01146526180179287, -0.007847949250506182, 0.1372692638829024, -1.240031967956483, 1.31846819961503, -0.2644632250616801, 0.09174663255667607, -0.08263104799555555, -0.06460783081261469, 0.004410425919287509, -0.007847949250506178, 0.005371905893059449, -0.09396054230867185, 0.8487994536684709, -0.9024889006343455, 0.1810245596471862, -0.06280039031442253, 0.05656079053370189, 0.04422393366747207, -0.07714320009289319, 0.1372692638829023, -0.0939605423086719, 1.643473226577989, -14.84643598858208, 15.78552346583742, -3.166318646351221, 1.098447897004905, -0.9893104342762569, -0.773525239112454, 0.6968787587238311, -1.240031967956483, 0.8487994536684713, -14.84643598858208, 134.1163689182897, -142.5996845533002, 28.60317183285527, -9.922909680513776, 8.937008402565429, 6.987697007968716, -0.7409587059910465, 1.31846819961503, -0.9024889006343462, 15.78552346583741, -142.5996845533002, 151.6195987015546, -30.41241955390257, 10.55056740429831, -9.502304523562367, -7.429692565696146, 0.1486242361257839, -0.2644632250616801, 0.1810245596471862, -3.166318646351221, 28.60317183285527, -30.41241955390257, 6.100235530521243, -2.116271809047857, 1.906007365633279, 1.490275197926278, -0.05156018640273603, 0.09174663255667606, -0.06280039031442254, 1.098447897004905, -9.922909680513776, 10.55056740429831, -2.116271809047857, 0.7341694181090099, -0.6612252322956818, -0.5170009212455956, 0.04643736907370834, -0.08263104799555554, 0.05656079053370192, -0.989310434276257, 8.937008402565427, -9.502304523562367, 1.906007365633279, -0.6612252322956818, 0.5955284938871694, 0.4656337431327625, 0.03630860018450286, -0.06460783081261467, 0.04422393366747208, -0.7735252391124541, 6.987697007968716, -7.429692565696146, 1.490275197926278, -0.5170009212455956, 0.4656337431327625, 0.3640712156837718},
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mortar_restrict_flatten(const int i, const int j, const int k, const int n, 
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
