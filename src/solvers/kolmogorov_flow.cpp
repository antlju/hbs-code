#include "includes.h"
#include <chrono>

using Clock=std::chrono::high_resolution_clock;

void set_force(Pencil &force, const SolverParams &params, const Grid &grid, const Int j)
{
        Real Amp = 1.0;
        Real kf = params.kf;
        for (size_t i=0;i<force.nx_;i++)
        {
                force(i,0,0) = Amp*sin(kf*grid.x(j));
                force(i,0,1) = 0.0;
                force(i,0,2) = 0.0;
        }
        
} //End set_force

/// Function to calculate the source term for Poisson eq: rho/(2*alpha^k*dt)*div(u*)
void calc_divustar(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        size_t Nx = meshCntr.u.nx_;
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        /// Parameters
        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho=params.rho;
        Real dt=params.dt;
        
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
                        //Copy from mesh to bundle
                        ff2bundle(meshCntr.ustar,ustarBndl,j,k);

                        // Calculate div(u*) and multiply by proper factor
                        // for use as source in Poisson eq.
                        div(ustarBndl,divustarPncl,xfac,yfac,zfac);
                        divustarPncl = divustarPncl*(rho/(2*alphak*dt));

                        // Copy from pencil to mesh
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
                                - gradpPncl*(2*alphak*dt/rho);

                        //Copy from pencil to mesh
                        pencil2ff(ustarPncl,meshCntr.ustar,j,k);
                        
                }
        }
}
void calc_pncl_diagnostics(const Pencil &u, const Pencil &omega, const Pencil &divu, Stats &stats)
{
        stats.calc_pncl_absmax(u);
	stats.calc_pncl_omega2(omega);
	stats.calc_pncl_Pavgs(divu);
        //stats.calc_pncl_rms(u);
        //std::cout << "umax: " << stats.umax << std::endl;
        //std::cout << "urms: " << stats.urms << std::endl;
}
/// Function to calculate RHSk = (u_j Del_j u_i)-nu*Lapl(u) at some RK substep k.
/// Here we also calculate necessary 
void calc_RHSk(MeshContainer &meshCntr, SolverParams &params, Stats &stats, const Grid &grid, const Int k_rk)
{
        size_t Nx = meshCntr.u.nx_;
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho = params.rho;
        
        /// Create pencils and bundles;
        Bundle uBndl(Nx,3);
        Pencil uPncl(Nx,3);
        Pencil uDeluPncl(Nx,3);
        Pencil LapluPncl(Nx,3);
        Pencil rhskPncl(Nx,3);
        Pencil forcePncl(Nx,3);
	Pencil omegaPncl(Nx,3); // omega = curl(u)
	Pencil divuPncl(Nx,1);
        
        //Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle
                        ff2bundle(meshCntr.u,uBndl,j,k);
                        ff2pencil(meshCntr.u,uPncl,j,k);

                        
                        //If k_rk = 1: Calculate pencil diagnostics from u
                        if (k_rk == 1)
			{
				/// Statistics set
				stats.umax_old = stats.umax;
				stats.urms_old = stats.urms;
				stats.umax = 0.0;
				stats.urms = 0.0;
				stats.omega2 = 0.0;
				stats.P = 0.0;
				stats.P2 = 0.0;
				
				//calculate omega = curl(u)
				div(uBndl,divuPncl,xfac,yfac,zfac);
				curl(uBndl,omegaPncl,xfac,yfac,zfac);
                                calc_pncl_diagnostics(uPncl,omegaPncl,divuPncl,stats);
			}
                        
                        
                        //Compute (u.grad)u on the bundle.
                        udotgradu(uBndl,uDeluPncl,xfac,yfac,zfac);

                        //Compute vector laplacian on the bundle.
                        vlapl(uBndl,LapluPncl,xfac*xfac,yfac*yfac,zfac*zfac);
                        LapluPncl = LapluPncl*(1.0/rho);

                        //Compute force pencil
                        set_force(forcePncl,params,grid,j);
                        
                        //Set RHSk pencil
                        rhskPncl = LapluPncl-uDeluPncl+forcePncl;

                        //Copy from pencil to mesh
                        pencil2ff(rhskPncl,meshCntr.RHSk,j,k);
			
                }
        }
	
        
        //std::cout << "urms final: " << stats.urms << std::endl;
        
}

/// Calls the Poisson solver
void PoissonEq(MeshContainer &meshCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        /// Since this is a inplace real-to-complex transform solver,
        /// here the input mesh psi acts as the input source as well as the result.
        PoissonSolve3D(meshCntr.psi,meshCntr.fftwPsi,0,grid.xlen);
        
}

/// Updates pressure p^k = p^(k-1) + psi^k
void update_pressure(MeshContainer &meshCntr)
{
        meshCntr.p = meshCntr.p + meshCntr.psi;
}

/// Calculate gradient of psi. Perhaps this should be done in Fourier space later.
void calc_gradpsi(MeshContainer &meshCntr, const SolverParams &params, const Grid &grid)
{
        size_t Nx = meshCntr.psi.nx_;
        size_t Ny = meshCntr.psi.ny_;
        size_t Nz = meshCntr.psi.nz_;
        
        Bundle psiBndl(Nx,1);
        Pencil gradpsiPncl(Nx,3);

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;

        // Loop over mesh (y,z)-plane 
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

/// Enforces solenoidal condition on velocity field by subtracting grad(psi) from ustar.
void enforce_solenoidal(MeshContainer &meshCntr, const SolverParams &params, const Grid &grid, const Int k_rk)
{
        ///Set param constants
        Real dt = params.dt, rho = params.rho;
        
        /// Create RK3 coefficient
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);

        apply_pbc(meshCntr.psi);
        calc_gradpsi(meshCntr,params,grid);
        meshCntr.u = meshCntr.ustar-(meshCntr.gradpsi*(2*alphak*dt/rho));
}

void update_timestep(SolverParams &params, Stats &stats, const Grid &grid)
{
        Real c1=0.5,c2=0.5; //Courant numbers for advection and diffusion respectively

        Real nu = params.viscosity,dx = grid.dx,L=grid.dx;
        
        //Factor of 1/3 since dx=dy=dz
        if (stats.umax < 1e-10)
                stats.umax = stats.umax_old;
        
        Real adv = (1.0/3)*c1*dx/stats.umax;
        Real diff = (1.0/3)*c2*pow(L,2)/nu;
        std::cout << "adv : " << adv << "\t diff: " << diff << std::endl;
        //Set new time step size according to CFL condition
        params.dt = std::min(adv,diff);

}

/// Function that performs RK3 substeps from k=1 to k=3
void RK3_stepping(MeshContainer &meshCntr, SolverParams &params, Stats &stats, const Grid &grid)
{
        
        /// From the previous step we have k=0 (time step n) data.
        /// We want to arrive at data for k=3 (time step n+1).
        /// (compare with Rosti & Brandt 2017)

        for (Int k_rk = 1;k_rk<=3;k_rk++)
        {
                /// First calculate RHSk = -D_j u_i u_j+(nu/rho)*Lapl(u))+force
                /// Then calc. u* = u+(2*dt*(alpha(k)/rho))*grad(p)
                ///                        +(dt*beta(k))*RHSk+(dt*gamma(k))*RHSk_1
                /// This is all done within the bundle/pencil framework.
                meshCntr.RHSk_1 = meshCntr.RHSk;

                apply_pbc(meshCntr.u);
                calc_RHSk(meshCntr,params,stats,grid,k_rk); //diagnostics also calculated here
                // If k_rk == 1 update the timestep dt
                if (k_rk == 1)
                        update_timestep(params,stats,grid);
                
                calc_ustar(meshCntr,params,grid,k_rk);

                /// Solve the Poisson eq for Psi.
                apply_pbc(meshCntr.ustar);
                calc_divustar(meshCntr,params,grid,k_rk);
                PoissonEq(meshCntr,params,grid,k_rk);

                /// Update pressure and enforce solenoidal condition
                update_pressure(meshCntr);
                enforce_solenoidal(meshCntr,params,grid,k_rk);
        }
        
} //End RK3_stepping()

void calc_Re(SolverParams &params, const Stats &stats, const Grid &grid)
{
        //std::cout << stats.urms << "\t" << grid.dx << "\t" << params.viscosity << std::endl;
        //params.Re = stats.urms*grid.dx/params.viscosity;
        //std::cout << "Re: " << params.Re << "\t sqrt(2): " << sqrt(2.0) << std::endl;
}

void calc_mesh_avg(const MeshContainer &meshCntr, Stats &stats, const Grid &grid)
{
	Real volsize = meshCntr.u.nx_*meshCntr.u.ny_*meshCntr.u.nz_;
	
	Real avgomega2 = stats.omega2/volsize;
	Real avgP = stats.P/volsize;
	Real avgP2 = stats.P2/volsize;
	//std::cout << "<Omega^2>: " << avgomega2 << std::endl;
	//std::cout << "<P>: " << avgP << std::endl;
	//std::cout << "<P^2>: " << avgP2 << std::endl;
	
}

Int main()
{
        auto t1 = Clock::now();
 
        /// Time settings.
        const Int maxtsteps = 100;

        /// Create parameter object and initialise parameters.
        SolverParams params;
        params.maxTimesteps = maxtsteps;
        params.currentTimestep = 0;
        params.kf = 1.0; //Kolmogorov frequency
        params.rho = 1.0;
        params.viscosity = 1.0;
        
        /// Set grid sizes
        const Real L0 = -M_PI, L1 = M_PI; // x,y,z in [-pi,pi]
        const Int Nsize = 256;
        
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

        /// At first we get initial timestep size from forcing
        Pencil initforcePncl(Nsize,1);
        for (size_t i=0;i<uu.nx_;i++)
                initforcePncl(i,0,0) = sin(params.kf*grid.x(i));

	/// Initial set of the timesteps
        stats.calc_pncl_absmax(initforcePncl);
        stats.umax_old = stats.umax;
        update_timestep(params,stats,grid);
        
        //std::cout << 0 << " : " << params.dt << std::endl;
        
/// Run Runge-Kutta 3 substepping
        for (Int ts = 1;ts<params.maxTimesteps;ts++)
        {
                RK3_stepping(meshCntr,params,stats,grid);
		//calc_mesh_avg(meshCntr,stats,grid);
                //stats.urms = stats.calc_mesh_rms(meshCntr.u);
                //calc_Re(params,stats,grid);
                //std::cout << ts << " : " << params.dt << std::endl;
        }

	
        // Write u_0(x,y,z) at last step to file
        writeToFile_1DArr(set_fname("kolm_rk3_x",".dat",params.maxTimesteps),
                          meshCntr.u,0,grid);

	writeStatsToFile(set_fname("kolm_rk3_x",".stats",params.maxTimesteps),meshCntr,stats,grid);
        
        auto t2 = Clock::now();
        std::cout << "Delta t2-t1: " 
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds" << std::endl;
        
        return 0;
        
} //End main()
