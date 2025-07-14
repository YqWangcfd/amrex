
#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <fstream>
#include <iomanip>

namespace amrex {

std::string LevelPath (int level, const std::string &levelPrefix)
{
    return Concatenate(levelPrefix, level, 1);  // e.g., Level_5
}

std::string MultiFabHeaderPath (int level,
                                const std::string &levelPrefix,
                                const std::string &mfPrefix)
{
    return LevelPath(level, levelPrefix) + '/' + mfPrefix;  // e.g., Level_4/Cell
}

std::string LevelFullPath (int level,
                           const std::string &plotfilename,
                           const std::string &levelPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
        r += '/';
    }
    r += LevelPath(level, levelPrefix);  // e.g., plt00005/Level_5
    return r;
}

std::string MultiFabFileFullPrefix (int level,
                                    const std::string& plotfilename,
                                    const std::string &levelPrefix,
                                    const std::string &mfPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
        r += '/';
    }
    r += MultiFabHeaderPath(level, levelPrefix, mfPrefix);
    return r;
}


void
PreBuildDirectorHierarchy (const std::string &dirName,
                           const std::string &/*subDirPrefix*/,
                           int nSubDirs, bool callBarrier)
{
  UtilCreateCleanDirectory(dirName, false);  // ---- dont call barrier
  for(int i(0); i < nSubDirs; ++i) {
    const std::string &fullpath = LevelFullPath(i, dirName);
    UtilCreateCleanDirectory(fullpath, false);  // ---- dont call barrier
  }

  if(callBarrier) {
    ParallelDescriptor::Barrier();
  }
}


void
WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                            int nlevels,
                            const Vector<BoxArray> &bArray,
                            const Vector<std::string> &varnames,
                            const Vector<Geometry> &geom,
                            Real time,
                            const Vector<int> &level_steps,
                            const Vector<IntVect> &ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix)
{
//        BL_PROFILE("WriteGenericPlotfileHeader()");

        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

        HeaderFile.precision(17);

        // ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (const auto & varname : varnames) {
            HeaderFile << varname << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << geom[i].Domain() << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
            }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
            HeaderFile << level_steps[level] << '\n';

            const IntVect& domain_lo = geom[level].Domain().smallEnd();
            for (int i = 0; i < bArray[level].size(); ++i)
            {
                // Need to shift because the RealBox ctor we call takes the
                // physical location of index (0,0,0).  This does not affect
                // the usual cases where the domain index starts with 0.
                const Box& b = amrex::shift(bArray[level][i], -domain_lo);
                // auto const& dx0 = geom[level].CellSize();
                // Real dx[AMREX_SPACEDIM];
                // Real gaussArray[5] = {0.5, 1.0, 2.0, 1.0, 0.5};
                // for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                //     dx[j] = dx0[j]*gaussArray[i%5];
                // }
                // RealBox loc = RealBox(b, dx, geom[level].ProbLo());
                RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
                }
            }

            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
        }
}


void
WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                         const Vector<const MultiFab*>& mf,
                         const Vector<std::string>& varnames,
                         const Vector<Geometry>& geom, Real time,
                         const Vector<int>& level_steps,
                         const Vector<IntVect>& ref_ratio,
                         const std::string &versionName,
                         const std::string &levelPrefix,
                         const std::string &mfPrefix,
                         const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            std::string ed = plotfilename;
            ed.append("/").append(d);
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) { FileOpenFailed(HeaderFileName); }
            WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                       geom, time, level_steps, ref_ratio, versionName,
                                       levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
        }
    }
}

// write a plotfile to disk given:
// -plotfile name
// -vector of MultiFabs
// -vector of Geometrys
// variable names are written as "Var0", "Var1", etc.
// refinement ratio is computed from the Geometry vector
// "time" and "level_steps" are set to zero
void WriteMLMF (const std::string &plotfilename,
                const Vector<const MultiFab*>& mf,
                const Vector<Geometry> &geom)
{
    int nlevs = static_cast<int>(mf.size());
    int ncomp = mf[0]->nComp();

    // variables names are "Var0", "Var1", etc.
    Vector<std::string> varnames(ncomp);
    for (int i=0; i<ncomp; ++i) {
        varnames[i] = "Var" + std::to_string(i);
    }

    // compute refinement ratio by looking at hi coordinates of domain at each level from
    // the geometry object
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int i = 0; i < nlevs-1; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int rr = (geom[i+1].Domain()).bigEnd(d)/(geom[i].Domain()).bigEnd(d);
            ref_ratio[i][d] = rr;
        }
    }

    // set step_array to zero
    Vector<int> step_array(nlevs,0);

    // set time to zero
    Real time = 0.;

    WriteMultiLevelPlotfile(plotfilename, nlevs, mf, varnames,
                            geom, time, step_array, ref_ratio);

}


void
WriteMultiLevelPlotfileHeaders (const std::string & plotfilename, int nlevels,
                                const Vector<const MultiFab*> & mf,
                                const Vector<std::string>     & varnames,
                                const Vector<Geometry>        & geom,
                                Real time, const Vector<int>  & level_steps,
                                const Vector<IntVect>     & ref_ratio,
                                const std::string         & versionName,
                                const std::string         & levelPrefix,
                                const std::string         & mfPrefix,
                                const Vector<std::string> & extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            std::string ed = plotfilename;
            ed.append("/").append(d);
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                   geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }

    for (int level = 0; level <= finest_level; ++level) {
        const MultiFab * data;
        data = mf[level];
        VisMF::WriteOnlyHeader(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

}

void
WriteSingleLevelPlotfile (const std::string& plotfilename,
                          const MultiFab& mf, const Vector<std::string>& varnames,
                          const Geometry& geom, Real time, int level_step,
                          const std::string &versionName,
                          const std::string &levelPrefix,
                          const std::string &mfPrefix,
                          const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}




#ifdef AMREX_USE_EB
void
EB_WriteSingleLevelPlotfile (const std::string& plotfilename,
                             const MultiFab& mf, const Vector<std::string>& varnames,
                             const Geometry& geom, Real time, int level_step,
                             const std::string &versionName,
                             const std::string &levelPrefix,
                             const std::string &mfPrefix,
                             const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    EB_WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                               level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
EB_WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                            const Vector<const MultiFab*>& mf,
                            const Vector<std::string>& varnames,
                            const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                            const Vector<IntVect>& ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix,
                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf[0]->hasEBFabFactory(),
                                     "EB_WriteMultiLevelPlotfile: does not have EB Factory");

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            std::string ed = plotfilename;
            ed.append("/").append(d);
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        Vector<std::string> vn = varnames;
        vn.push_back("vfrac");
        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, vn,
                                   geom, time, level_steps, ref_ratio, versionName,
                                   levelPrefix, mfPrefix);

        for (int lev = 0; lev < nlevels; ++lev) {
            HeaderFile << "1.0e-6\n";
        }
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const int nc = mf[level]->nComp();
        MultiFab mf_tmp(mf[level]->boxArray(),
                        mf[level]->DistributionMap(),
                        nc+1, 0);
        MultiFab::Copy(mf_tmp, *mf[level], 0, 0, nc, 0);
        auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf[level]->Factory());
        MultiFab::Copy(mf_tmp, factory.getVolFrac(), 0, nc, 1, 0);
        VisMF::Write(mf_tmp, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

#endif
}


namespace amrex
{
void WriteMultiLevelPlotfile_SD(const std::string& plotfilename, int nlevels,
                                const Vector<const MultiFab*>& mf,
                                const Vector<std::string>& varnames,
                                const Vector<Geometry>& geom, Real time,
                                const Vector<int>& level_steps,
                                const Vector<IntVect>& ref_ratio,
                                const std::string &versionName,
                                const std::string &levelPrefix,
                                const std::string &mfPrefix,
                                const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile_SD()");
    
    int finest_level = nlevels - 1;

    constexpr int ORDER = 5;

    // Vector of MultiFabs for mapping cell-based data
    Vector<MultiFab> mf_sd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray ba_sd = mf[lev]->boxArray();
        int ncomps = mf[lev]->nComp()/std::pow(ORDER,AMREX_SPACEDIM);
        // Print() << "ncomps of mf_sd" << ncomps << std::endl;
        ba_sd.refine(ORDER);
        // Print() << "ba_sd= " << ba_sd << std::endl;
        mf_sd[lev].define(ba_sd, mf[lev]->DistributionMap(), ncomps, 0);
    }

    for (int lev(0); lev <= finest_level; ++lev) {
        MultiFab const& S_old = *mf[lev];
        int ncomps = S_old.nComp()/std::pow(ORDER,AMREX_SPACEDIM);
        for (MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& fab = S_old.const_array(mfi);
            Array4<Real> const& fab_sd = mf_sd[lev].array(mfi);

            ParallelFor(bx, ncomps,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
            {
                int zloop = (AMREX_SPACEDIM==2)? 1:ORDER;
                for (int iz(0); iz < zloop; iz++) {
                    for (int iy(0); iy < ORDER; iy++) {
                        for (int ix(0); ix < ORDER; ix++) {
                            int ii = i*ORDER + ix;
                            int jj = j*ORDER + iy;
                            int kk = (AMREX_SPACEDIM == 2)? 
                            int(0) : k*ORDER + iz;

                            // Mapping from coarse box with refined components
                            // to refined box with coarse components
                            int pid = (AMREX_SPACEDIM == 2)?
                            iy*ORDER+ix : iz*ORDER*ORDER+iy*ORDER+ix;
                            int nid = n * std::pow(ORDER, AMREX_SPACEDIM);
                            fab_sd(ii,jj,kk,n) = fab(i,j,k,nid+pid);
                            // Print() << "fab(0,0,0,9*25)= " <<
                            // fab(0,0,0,225) << std::endl;
                            // Print() << "fab= " << fab(i,j,k,0) << std::endl;
                        }
                    }
                }
            });
        }
    }

    // build new geometry for mesh refinement
    Vector<Geometry> geom_sd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        // 获取原始几何参数
        const auto& dx = geom[lev].CellSizeArray();
        const auto& prob_lo = geom[lev].ProbLoArray();
        const auto& prob_hi = geom[lev].ProbHiArray();
        
        // 创建新的几何参数，网格间距缩小5倍
        Array<Real,AMREX_SPACEDIM> dx_sd;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            dx_sd[i] = dx[i] / ORDER;  // ORDER = 5
        }
        
        // 定义新的box
        Box domain = geom[lev].Domain();
        domain.refine(ORDER);
        
        // 创建新的geometry
        geom_sd[lev].define(domain);
    }

    // Vector of MultiFabs for nodal data
    Vector<MultiFab> mf_nd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray nodal_grids = mf_sd[lev].boxArray(); 
        nodal_grids.surroundingNodes();
        mf_nd[lev].define(nodal_grids, mf_sd[lev].DistributionMap(), 3, 0);
        mf_nd[lev].setVal(0.);
    }

    // Fill flux points physical coordinates MF
    Real x_Legendre[6] = {-1.0,
                      -0.8611363115940526,
                      -0.33998104358485626,
                       0.33998104358485626,
                       0.8611363115940526};
    for (int lev(0); lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();
        // Print() << "ProbloArray[0]" << prob_lo[0] << std::endl;
        for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            // Print() << "mfi.tilebox()= " << mfi.tilebox() << std::endl;
            Array4<Real> mf_arr = mf_nd[lev].array(mfi);
            ParallelFor(bx, 
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                int ii = i%5, jj = j%5, kk = k%5;
                int ix = i/5, iy = j/5, iz = k/5;
                Real x = prob_lo[0] + (ix+0.5)*dx[0];
                Real y = prob_lo[1] + (iy+0.5)*dx[1];
                Real z = (AMREX_SPACEDIM==2)? 0.0:
                prob_lo[2] + (iz+0.5)*dx[2];

                mf_arr(i,j,k,0) = dx[0]*x_Legendre[ii]/2.0 + x;
                mf_arr(i,j,k,1) = dx[1]*x_Legendre[jj]/2.0 + y;
                mf_arr(i,j,k,2) = (AMREX_SPACEDIM==2)? 0.0:
                dx[2]*x_Legendre[kk]/2.0 + z;
                
                // convert from coordinates to displacement
                mf_arr(i,j,k,0) -= prob_lo[0] + i*dx[0]/ORDER;
                mf_arr(i,j,k,1) -= prob_lo[1] + j*dx[1]/ORDER;
                mf_arr(i,j,k,2) -= (AMREX_SPACEDIM==2)? 0.0:
                prob_lo[2] + k*dx[2]/ORDER;

            });
        }
    }


    WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf_sd),
                                        GetVecOfConstPtrs(mf_nd),
                                        varnames,
                                        geom_sd, time, level_steps, ref_ratio);

}

void
WriteSingleLevelPlotfile_SD (const std::string& plotfilename,
                            const MultiFab& mf, const Vector<std::string>& varnames,
                            const Geometry& geom, Real time, int level_step,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix,
                            const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfile_SD(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}


void
WriteMultiLevelPlotfileWithTerrain (const std::string& plotfilename, int nlevels,
                                    const Vector<const MultiFab*>& mf,
                                    const Vector<const MultiFab*>& mf_nd,
                                    const Vector<std::string>& varnames,
                                    const Vector<Geometry>& my_geom,
                                    Real time,
                                    const Vector<int>& level_steps,
                                    const Vector<IntVect>& rr,
                                    const std::string &versionName,
                                    const std::string &levelPrefix,
                                    const std::string &mfPrefix,
                                    const Vector<std::string>& extra_dirs) 
{
    BL_PROFILE("WriteMultiLevelPlotfileWithTerrain()");

    int finest_level = nlevels - 1;

    AMREX_ALWAYS_ASSERT(nlevels <= mf.size());
    AMREX_ALWAYS_ASSERT(nlevels <= rr.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());
    // Print() << "mf[0]->nComp() = " << mf[0]->nComp() << ", varnames.size() = "
    // << varnames.size() << std::endl;
    AMREX_ALWAYS_ASSERT(mf[0]->nComp() == varnames.size());

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeaderWithTerrain(HeaderFile, nlevels, boxArrays, varnames,
                                                  my_geom, time, level_steps, rr, versionName,
                                                  levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "Nu_nd";
    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data        , MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
            VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
        }
    }
}

void
WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                            int nlevels,
                                            const Vector<BoxArray> &bArray,
                                            const Vector<std::string> &varnames,
                                            const Vector<Geometry>& my_geom,
                                            Real my_time,
                                            const Vector<int>& level_steps,
                                            const Vector<IntVect>& my_ref_ratio,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix)
{
    AMREX_ALWAYS_ASSERT(nlevels <= bArray.size());
    AMREX_ALWAYS_ASSERT(nlevels <= my_ref_ratio.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());

    int finest_level  = nlevels - 1;

    HeaderFile.precision(17);

    // ---- this is the generic plot file type name
    HeaderFile << versionName << '\n';

    HeaderFile << varnames.size() << '\n';

    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        HeaderFile << varnames[ivar] << "\n";
    }
    HeaderFile << AMREX_SPACEDIM << '\n';
    HeaderFile << my_time << '\n';
    HeaderFile << finest_level << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbLo(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbHi(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < finest_level; ++i) {
        HeaderFile << my_ref_ratio[i][0] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << my_geom[i].Domain() << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << level_steps[i] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            HeaderFile << my_geom[i].CellSize()[k] << ' ';
        }
        HeaderFile << '\n';
    }
    HeaderFile << (int) my_geom[0].Coord() << '\n';
    HeaderFile << "0\n";

    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << level << ' ' << bArray[level].size() << ' ' << my_time << '\n';
        HeaderFile << level_steps[level] << '\n';

        const IntVect& domain_lo = my_geom[level].Domain().smallEnd();
        for (int i = 0; i < bArray[level].size(); ++i)
        {
            // Need to shift because the RealBox ctor we call takes the
            // physical location of index (0,0,0).  This does not affect
            // the usual cases where the domain index starts with 0.
            const Box& b = shift(bArray[level][i], -domain_lo);
            RealBox loc = RealBox(b, my_geom[level].CellSize(), my_geom[level].ProbLo());
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
        }

        HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
    }
    HeaderFile << "1" << "\n";
    HeaderFile << "3" << "\n";
    HeaderFile << "amrexvec_nu_x" << "\n";
    HeaderFile << "amrexvec_nu_y" << "\n";
    HeaderFile << "amrexvec_nu_z" << "\n";
    std::string mf_nodal_prefix = "Nu_nd";
    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
    }
}
}
