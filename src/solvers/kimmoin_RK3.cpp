#include "includes.h"
void calc_divustar(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        size_t Nx = meshCntr.u.nx_;
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        /// Parameters
        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho = params.rho;
        Real dt = params.dt;
        
        /// Create RK3 coefficients
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);

        /// Create pencils and bundles;
        Bundle ustarBndl(Nx,3);
        Pencil divustarPncl(Nx,1);

        //Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        ff2bundle(meshCntr.ustar,ustarBndl,j,k);
                        
                        div(ustarBndl,divustarPncl,xfac,yfac,zfac);
                        divustarPncl = divustarPncl*(rho/(2*alphak*dt));
                        
                        pencil2ff(divustarPncl,meshCntr.psi,j,k);
                        
                }
        }
}

/// Calculate u* = u+(2*dt*(alpha(k)/rho))*grad(p)+(dt*beta(k))*RHSk+(dt*gamma(k))*RHSk_1
void calc_ustar(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        size_t Nx = meshCntr.u.nx_;
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        /// Parameters
        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho = params.rho;
        Real dt = params.dt;
        
        /// Create RK3 coefficients
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);
        Real betak = coeffs.beta(k_rk);
        Real gammak = coeffs.beta(k_rk);

        /// Create pencils and bundles;
        Bundle pBndl(Nx,1);
        Pencil uPncl(Nx,3);
        Pencil gradpPncl(Nx,3);
        Pencil ustarPncl(Nx,3);
        Pencil rhskPncl(Nx,3);
        Pencil rhsk_1Pncl(Nx,3);

        //Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle/pencil
                        ff2bundle(meshCntr.p,pBndl,j,k);
                        ff2pencil(meshCntr.u,uPncl,j,k);
                        ff2pencil(meshCntr.RHSk,rhskPncl,j,k);
                        ff2pencil(meshCntr.RHSk_1,rhsk_1Pncl,j,k);
                        
                        //Compute gradient of the pressure (scalar bundle -> vector pencil)
                        sgrad(pBndl,gradpPncl,xfac,yfac,zfac);

                        // Set ustar
                        ustarPncl = uPncl + rhskPncl*(dt*betak) + rhsk_1Pncl*(dt*gammak)
                                + gradpPncl*(2*alphak*dt/rho);

                        //Copy from pencil to mesh
                        pencil2ff(ustarPncl,meshCntr.ustar,j,k);
                        
                }
        }
}

/// Function to calculate RHSk = (u_j Del_j u_i)-nu*Lapl(u) at some RK substep k.
void calc_RHSk(MeshContainer &meshCntr, SolverParams &params, const Grid &grid)
{
        size_t Nx = meshCntr.u.nx_;
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho = params.rho;
        
        /// Create pencils and bundles;
        Bundle uBndl(Nx,3);
        Pencil uDelu(Nx,3);
        Pencil Laplu(Nx,3);
        Pencil rhskPncl(Nx,3);

        //Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle
                        ff2bundle(meshCntr.u,uBndl,j,k);
                        
                        //Compute (u.grad)u on the bundle.
                        udotgradu(uBndl,uDelu,xfac,yfac,zfac);

                        //Compute vector laplacian on the bundle.
                        vlapl(uBndl,Laplu,xfac*xfac,yfac*yfac,zfac*zfac);
                        Laplu = Laplu*(1.0/rho);

                        //Set RHSk pencil
                        rhskPncl = uDelu-Laplu;

                        //Copy from pencil to mesh
                        pencil2ff(rhskPncl,meshCntr.RHSk,j,k);
                                        
                }
        }
        
}

void PoissonEq(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        /// Since this is a inplace real-to-complex transform solver,
        /// here the input mesh psi acts as the input source as well as the result.
        PoissonSolve3D(meshCntr.psi,meshCntr.fftwPsi,0,grid.xlen);
        
}

void update_pressure(MeshContainer &meshCntr)
{
        meshCntr.p = meshCntr.p + meshCntr.psi;
}

void calc_gradpsi(MeshContainer &meshCntr, const SolverParams &params, const Grid &grid)
{
        size_t Nx = meshCntr.psi.nx_;
        size_t Ny = meshCntr.psi.ny_;
        size_t Nz = meshCntr.psi.nz_;
        
        Bundle psiBndl(Nx,1);
        Pencil gradpsiPncl(Nx,3);

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        ff2bundle(meshCntr.psi,psiBndl,j,k); //Copy
                        sgrad(psiBndl,gradpsiPncl,xfac,yfac,zfac); //Compute gradient on bundle
                        pencil2ff(gradpsiPncl,meshCntr.gradpsi,j,k);
                        
                }
        }
}

void enforce_solenoidal(MeshContainer &meshCntr, const SolverParams &params, const Grid &grid, const Int k_rk)
{
        ///Set param constants
        Real dt = params.dt, rho = params.rho;
        
        /// Create RK3 coefficient
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);

        calc_gradpsi(meshCntr,params,grid);
        meshCntr.u = meshCntr.ustar-(meshCntr.gradpsi*(2*alphak*dt/rho));
}

void RK3_timestepping(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, Stats &stats)
{
        /// From the previous step we have k=0 (time step n) data.
        /// We want to arrive at data for k=3 (time step n+1).

        /// Let k=1.
        Int k_rk = 1;
        /// First calculate RHSk = D_j u_i u_j-(nu/rho)*Lapl(u) )
        /// Then calc. u* = u+(2*dt*(alpha(k)/rho))*grad(p)
        ///                        +(dt*beta(k))*RHSk+(dt*gamma(k))*RHSk_1
        /// This is all done within the bundle/pencil framework.
        meshCntr.RHSk_1 = meshCntr.RHSk;
        calc_RHSk(meshCntr,params,grid);
        calc_ustar(meshCntr,params,grid,k_rk);

        /// Solve the Poisson eq for Psi.
        calc_divustar(meshCntr,params,grid,k_rk);
        PoissonEq(meshCntr,params,grid,k_rk);

        /// Update pressure and enforce solenoidal condition
        update_pressure(meshCntr);
        enforce_solenoidal(meshCntr,params,grid,k_rk);
        
} //End RK3_timestepping()


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
        params.rho = 1.0;
        
        /// Set grid sizes
        const Real L0 = -M_PI, L1 = M_PI;
        const Int Nsize = 4;
        
        /// Create and initialise uniform 3D finite difference grid object.
        Grid grid(Nsize,Nsize,Nsize,L0,L1);

        /// Create run statistics object
        Stats stats;

        /// Create mesh objects
        Mesh uu(Nsize,Nsize,Nsize,3); /// velocity vector field
        Mesh ustar(Nsize,Nsize,Nsize,3); /// ustar
        Mesh RHSk(Nsize,Nsize,Nsize,3); /// RHS^k vector field
        Mesh RHSk_1(Nsize,Nsize,Nsize,3); /// RHS^(k-1) vector field
        Mesh pp(Nsize,Nsize,Nsize,1); /// pressure scalar field
        Mesh psi(Nsize,Nsize,Nsize,1); /// Psi scalar field
        Mesh gradpsi(Nsize,Nsize,Nsize,3); /// grad(psi) vector field
        fftwMesh fftwPsi(Nsize,Nsize,Nsize); /// FFTW3 memory scalar field
        
        /// Create container object with references to meshes
        MeshContainer meshCntr(uu,ustar,RHSk,RHSk_1,pp,psi,gradpsi,fftwPsi);

        RK3_timestepping(meshCntr,params,grid,stats);
        return 0;
        
}//End main()
