#include "includes.h"

//void init_velocity(Mesh &u, const Pencil &x, SolverParams &params)
void init_velocity(Mesh &u)
{
        /// Kolmogorov frequency
        //Real kf = params.kf;

        /// Set initial state
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                //Real term = sin(x(i))+sin(x(j))+sin(x(k));
                                //Real term = sin(kf*x(j));
                                u(i,j,k,0) = 0.0;
                                u(i,j,k,1) = 0.0;
                                u(i,j,k,2) = 0.0;
                        }
                }
        }
        
}//End init_velocity()


void set_force(Pencil &force, const Grid &grid, const SolverParams &params, const Int j)
{
        Real Amp = 1.0, kf = params.kf;

        for (size_t i=0;i<force.nx_;i++)
        {
                force(i,0,0) = Amp*sin(kf*grid.x(j));
                force(i,0,1) = 0.0;
                force(i,0,2) = 0.0;
        }
        
}

/// Computes u* = u_old+dt*(-grad(p)-(u dot grad)u+nu*Lapl(u)+force)
void compute_ustar(MeshContainer &meshcntr, SolverParams &params, Grid &grid, Stats &stats)
{
        // Sizes
        Int nx = meshcntr.u.nx_, ny = meshcntr.u.ny_, nz = meshcntr.u.nz_;
        
        /// Create bundle and pencil objects
        Bundle uBndl(nx,3); /// Vector bundle
        Bundle pBndl(nx,1); /// Scalar bundle

        Pencil ustar(nx,3); /// Vector pencil
        Pencil uDelu(nx,3); /// Vector pencil
        Pencil uPncl(nx,3); /// Vector pencil, to save "old u"

        Pencil force(nx,3); /// Vector pencil

        
        for (size_t j=0;j<ny;j++)
        {
                for (size_t k=0;k<nz;k++)
                {
                        /// Copy from mesh operations
                        ff2bundle(meshcntr.u,uBndl,j,k);
                        ff2bundle(meshcntr.p,pBndl,j,k);

                        // Copy "old" u to pencil
                        bundle2pencil(uBndl,uPncl);

                        set_force(force,grid,params,j);
    
                        ustar = uPncl + force;

                }
        }
        

}

Int main()
{

        /// Time settings.
        Real dt = 0.01;
        const Int maxtsteps = 10;

        /// Create parameter object and initialise parameters.
        SolverParams params;
        params.dt = dt;
        params.maxTimesteps = maxtsteps;
        params.currentTimestep = 0;
        params.kf = 0.1; //Kolmogorov frequency
        
        /// Set grid sizes
        const Real L0 = -M_PI, L1 = M_PI;
        const Int Nsize = 4;
        
        /// Initialise uniform 3D finite difference grid object.
        Grid grid(Nsize,Nsize,Nsize,L0,L1);

        /// Create run statistics object
        Stats stats;

        /// Create mesh objects
        Mesh u(Nsize,Nsize,Nsize,3); /// vector field
        Mesh du(Nsize,Nsize,Nsize,3); ///  vector field
        Mesh pp(Nsize,Nsize,Nsize,1); /// scalar field
        
        /// Create container object with references to meshes
        MeshContainer meshCntr(u,du,pp);

        compute_ustar(meshCntr,params,grid,stats);
        return 0;
        
}//End main()
