#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        std::vector<int> n_cell;
        std::vector<Real> prob_lo;
        std::vector<Real> prob_hi;
        int max_grid_size = 32;
        int n_iters=20;
        int n_out=1;
        amrex::Real B=-1.0;

        Vector<int> is_periodic(AMREX_SPACEDIM,1);

        ParmParse pp;
        pp.getarr("n_cell", n_cell);
        pp.getarr("prob_lo", prob_lo);
        pp.getarr("prob_hi", prob_hi);
        pp.query("max_grid_size", max_grid_size);
        pp.query("n_iters",n_iters);
        pp.query("n_out",n_out);

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        RealBox rb({prob_lo[0],prob_lo[1]}, 
                {prob_hi[0],prob_hi[1]});

        Box domain(IntVect{0,0},
                IntVect{n_cell[0]-1,n_cell[1]-1});

        geom.define(domain, &rb, CoordSys::cartesian, is_periodic.data());
        grids.define(domain); 
        grids.maxSize(max_grid_size); 
        dmap.define(grids); 

        MultiFab phi(grids, dmap, 1, 2);
        Array<MultiFab,AMREX_SPACEDIM> enrg;
        for(int d=0;d<AMREX_SPACEDIM;d++)
        {
           BoxArray ba = grids;
           ba.surroundingNodes(d);
           enrg[d].define(ba,dmap,1,0);
           enrg[d].setVal(0.0);
        }
        
        Array<MultiFab,AMREX_SPACEDIM> totalenrg;
        for(int d=0;d<AMREX_SPACEDIM;d++)
        {
           BoxArray ba = grids;
           ba.surroundingNodes(d);
           totalenrg[d].define(ba,dmap,1,0);
           totalenrg[d].setVal(0.0);
        }

        phi.setVal(1.0);
        amrex::Real vacfrac=0.2;
        for(MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> phi_arr = phi.array(mfi);

            amrex::ParallelForRNG(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine)
            {
                //random pattern
                /*if(amrex::Random(engine) < vacfrac)
                {
                    phi_arr(i,j,k)=0.0;
                }
                else
                {
                    phi_arr(i,j,k)=1.0;
                }*/

                //chessboard
                phi_arr(i,j,k)=(i+j)%2;
            });
        }
        
        MultiFab plotfile_mf(grids, dmap, 1, 0);
        const Real strt_evolve = amrex::second();
        for(int iter=0;iter<n_iters;iter++)
        {
            if(iter%n_out==0)
            {
                MultiFab::Copy(plotfile_mf, phi,0,0,1,0);
                const std::string& plotfilename = amrex::Concatenate("plt", iter, 5);
                WriteSingleLevelPlotfile(plotfilename, plotfile_mf, 
                        {"vacancies"}, geom, 0.0, 0);
            }
            
            phi.FillBoundary(geom.periodicity());
            for(int d=0;d<AMREX_SPACEDIM;d++)
            {
                enrg[d].setVal(0.0);
            }
            
            //find total energy
            for(MFIter mfi(phi); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> phi_arr  = phi.array(mfi);
                Box bx_x = convert(bx, IntVect::TheDimensionVector(0));
                Box bx_y = convert(bx, IntVect::TheDimensionVector(1));
                GpuArray<Array4<Real>, AMREX_SPACEDIM> tenrg_arr{totalenrg[0].array(mfi),
                    totalenrg[1].array(mfi)};

                amrex::ParallelFor(bx_x,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect ivf(i,j);
                    IntVect ivR(i,j);
                    IntVect ivL(i,j);
                    ivL[0]-=1;
                    tenrg_arr[0](ivf)=(int(phi_arr(ivL)) & int(phi_arr(ivR)))*B;
                });
                
                amrex::ParallelFor(bx_y,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect ivf(i,j);
                    IntVect ivR(i,j);
                    IntVect ivL(i,j);
                    ivL[1]-=1;
                    tenrg_arr[1](ivf)=(int(phi_arr(ivL)) & int(phi_arr(ivR)))*B;
                });
             }

             
            amrex::Print()<<"iter, before phi sum:\t"<<iter<<"\t"<<phi.sum(0,false)<<"\n";
            amrex::Print()<<"iter, energy:\t"<<iter<<"\t"<<
            totalenrg[0].sum(0,false)+totalenrg[1].sum(0,false)<<"\n";

            for(MFIter mfi(phi); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx,1);
            
                GpuArray<Array4<Real>, AMREX_SPACEDIM> enrg_arr{enrg[0].array(mfi),enrg[1].array(mfi)};
                Array4<Real> phi_arr  = phi.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect iv(i,j);

                    if(int(phi_arr(iv))==0)
                    {
                        //amrex::Print()<<"phi_arr:"<<phi_arr(iv)<<"\t"<<iv<<"\n";
                        //directions x or y
                        for(int mvdir=0;mvdir<2;mvdir++)
                        {
                            //left/right or bottom/top
                            for(int mv=0;mv<2;mv++)
                            {
                                int mvnum=2*mv-1;

                                IntVect iv_n(i,j);
                                iv_n[mvdir]+=mvnum;
                                if(int(phi_arr(iv_n))==1)
                                {
                                    amrex::Real enrg1=0.0;
                                    amrex::Real enrg2=0.0;

                                    enrg1 += (int(phi_arr(iv_n)) & int(phi_arr(iv_n[0]-1,iv_n[1],k)))*B;
                                    enrg1 += (int(phi_arr(iv_n)) & int(phi_arr(iv_n[0]+1,iv_n[1],k)))*B;
                                    enrg1 += (int(phi_arr(iv_n)) & int(phi_arr(iv_n[0],iv_n[1]-1,k)))*B;
                                    enrg1 += (int(phi_arr(iv_n)) & int(phi_arr(iv_n[0],iv_n[1]+1,k)))*B;

                                    //amrex::Print()<<"enrg1, phi_arr(iv_n):"<<enrg1<<"\t"<<phi_arr(iv_n)<<"\n";

                                    for(int d=0;d<2;d++)
                                    {
                                        for(int m=0;m<2;m++)
                                        {
                                            IntVect iv_enrg2(i,j);
                                            iv_enrg2[d]+=2*m-1;
                                            if(iv_enrg2[0]!=iv_n[0] || iv_enrg2[1]!=iv_n[1])
                                            {
                                                enrg2 += (1 & int(phi_arr(iv_enrg2)))*B;
                                            }
                                        }
                                    }
                                    //amrex::Print()<<"enrg1,enrg2:"<<enrg1<<"\t"<<enrg2<<"\n";

                                    IntVect iv_f(i,j);
                                    iv_f[mvdir]+=mv;
                                    enrg_arr[mvdir](iv_f)=enrg2-enrg1;
                                    //amrex::Print()<<"enrg_arr:"<<enrg_arr(i,j,k,2*mvdir+mv)<<"\t"<<enrg2<<"\t"<<enrg1<<"\n";
                                }
                            }
                        }
                    }

                });
            }

            Real minmove_enrg=1e20;
            IntVect minmove_ind1(0,0);
            IntVect minmove_ind2(0,0);
            //GpuArray<IntVect, AMREX_SPACEDIM> minmove_ind1;
            int minmove_dim=0;

            for(int idim=0;idim<AMREX_SPACEDIM;idim++)
            {
                Real minval=enrg[idim].min(0,0,false);
                IntVect minind=enrg[idim].minIndex(0,0);

                if(minval < minmove_enrg)
                {
                    minmove_dim=idim;
                    minmove_enrg=minval;
                    minmove_ind1=minind;
                    minmove_ind2=minind;
                }
            }
            
            //this is a hack to make sure 
            //we check both sides when at the 
            //physical boundary so that we dont 
            //loose a vacancy into the ghost layers
            if(minmove_dim==0)
            {
                int idim=minmove_dim;
                if(minmove_ind1[idim]==0)
                {
                   minmove_ind2[idim]=(minmove_ind1[idim]+n_cell[idim]);
                }
                if(minmove_ind1[idim]==n_cell[idim])
                {
                   minmove_ind2[idim]=(minmove_ind1[idim]-n_cell[idim]);
                }
            }
            
            if(minmove_dim==1)
            {
                int idim=minmove_dim;
                if(minmove_ind1[idim]==0)
                {
                   minmove_ind2[idim]=(minmove_ind1[idim]+n_cell[idim]);
                }
                if(minmove_ind1[idim]==n_cell[idim])
                {
                   minmove_ind2[idim]=(minmove_ind1[idim]-n_cell[idim]);
                }
            }

            if(minmove_enrg < 0.0)
            {
                for(MFIter mfi(phi); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> phi_arr = phi.array(mfi);
                    Box bx_f = convert(bx, IntVect::TheDimensionVector(minmove_dim));

                    amrex::ParallelFor(bx_f,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if( (i==minmove_ind1[0] && j==minmove_ind1[1]) || 
                        (i==minmove_ind2[0] && j==minmove_ind2[1]) )

                        {
                            IntVect iv(i,j);
                            IntVect iv_swap(i,j);
                            iv_swap[minmove_dim]-=1;

                            //amrex::Print()<<"bx,iv,iv_swap:"<<bx<<"\t"<<iv<<"\t"<<iv_swap<<"\n";
                            
                            Real tmp=phi_arr(iv);
                            phi_arr(iv)=phi_arr(iv_swap);
                            phi_arr(iv_swap)=tmp;
                        }
                    });
                }
                /*amrex::Print()<<"iter, min enrg,minmove_ind,min dir:\t"<<iter<<"\t"<<
                minmove_enrg<<"\t"<<
                minmove_ind1<<"\t"<<minmove_ind2<<"\t"<<minmove_dim<<"\n";*/
                amrex::Print()<<"iter, after phi sum:\t"<<iter<<"\t"<<phi.sum(0,false)<<"\n";
                amrex::Print()<<"===============\n";
            }
            else
            {
                amrex::Print()<<"achieved stable state\n";
                break;
            }

        }
        Real end_evolve = amrex::second() - strt_evolve;

        amrex::Print()<<"Total time:"<<end_evolve<<"\n";

        MultiFab::Copy(plotfile_mf, phi,0,0,1,0);
        WriteSingleLevelPlotfile("finalplt", plotfile_mf, 
                                 {"vacancies"}, geom, 0.0, 0);
    }

    amrex::Finalize();
}
