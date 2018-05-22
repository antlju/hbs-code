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
void calc_divustar(MeshContainer &meshCntr, PBContainer &pbCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        /// Parameters
        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho=params.rho;
        Real dt=params.dt;
        
        /// Create RK3 coefficients
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);

        //Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle
                        ff2bundle(meshCntr.ustar,pbCntr.ustarBndl,j,k);

                        // Calculate div(u*) and multiply by proper factor
                        // for use as source in Poisson eq.
                        div(pbCntr.ustarBndl,pbCntr.dsPncl,xfac,yfac,zfac);
                        pbCntr.dsPncl = pbCntr.dsPncl*(rho/(2*alphak*dt));

                        // Copy from pencil to mesh
                        pencil2ff(pbCntr.dsPncl,meshCntr.psi,j,k);
                }
        }
}

/// Calculate u* = u+(2*dt*(alpha(k)/rho))*grad(p)+(dt*beta(k))*RHSk+(dt*gamma(k))*RHSk_1
void calc_ustar(MeshContainer &meshCntr, PBContainer &pbCntr, SolverParams &params, const Grid &grid, const Int k_rk)
{
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

	//Loop over mesh (y,z)-plane
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle/pencil
                        ff2bundle(meshCntr.p,pbCntr.pBndl,j,k); //pressure field
                        ff2pencil(meshCntr.u,pbCntr.uPncl,j,k); //velocity field
                        ff2pencil(meshCntr.RHSk,pbCntr.rhskPncl,j,k); // Right hand sides
                        ff2pencil(meshCntr.RHSk_1,pbCntr.rhsk_1Pncl,j,k);
                        
                        //Compute gradient of the pressure (scalar bundle -> vector pencil)
			pbCntr.dvPncl.fillPencil(0.0);
                        sgrad(pbCntr.pBndl,pbCntr.dvPncl,xfac,yfac,zfac);

                        // Set ustar
			// pbCntr.ustarPncl = pbCntr.uPncl + pbCntr.rhskPncl*(dt*betak) + pbCntr.rhsk_1Pncl*(dt*gammak) - pbCntr.dvPncl*(2*alphak*dt/rho);

			
			pbCntr.ustarPncl += pbCntr.uPncl;
			pbCntr.ustarPncl += pbCntr.rhskPncl*(dt*betak);
			pbCntr.ustarPncl += pbCntr.rhsk_1Pncl*(dt*gammak);
			pbCntr.ustarPncl += pbCntr.dvPncl*(-2*alphak*dt/rho);
			
				

                        //Copy from pencil to mesh
                        pencil2ff(pbCntr.ustarPncl,meshCntr.ustar,j,k);
                        
                }
        }
}
void calc_pncl_diagnostics(const Pencil &u, const Pencil &omega, Stats &stats)
{
        stats.calc_pncl_absmax(u);
	stats.calc_pncl_omega2(omega);
	//stats.calc_pncl_Pavgs(divu);
        //stats.calc_pncl_rms(u);
        //std::cout << "umax: " << stats.umax << std::endl;
        //std::cout << "urms: " << stats.urms << std::endl;
}

/// Function to calculate RHSk = (u_j Del_j u_i)-nu*Lapl(u) at some RK substep k.
/// Here we also calculate necessary 
void calc_RHSk(MeshContainer &meshCntr, PBContainer &pbCntr, SolverParams &params, Stats &stats, const Grid &grid, const Int k_rk)
{
        size_t Ny = meshCntr.u.ny_;
        size_t Nz = meshCntr.u.nz_;

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;
        Real rho = params.rho;

	//std::cout << "umax before calc: " << stats.umax << std::endl;
        //Loop over mesh (y,z)-plane

        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        //Copy from mesh to bundle
                        ff2bundle(meshCntr.u,pbCntr.uBndl,j,k);
                        ff2pencil(meshCntr.u,pbCntr.uPncl,j,k);

                        
                        //If k_rk = 1: Calculate pencil diagnostics from u
                        if (k_rk == 1)
			{
				
				/// Statistics set
				stats.umax_old = stats.umax;
				stats.urms_old = stats.urms;
				//stats.umax = 0.0;
				//stats.urms = 0.0;
				stats.omega2 = 0.0;
				stats.P = 0.0;
				stats.P2 = 0.0;
				
				//calculate omega = curl(u)
				
				//div(uBndl,divuPncl,xfac,yfac,zfac);
				//curl(pbCntr.uBndl,pbCntr.dvPncl,xfac,yfac,zfac);
                                //calc_pncl_diagnostics(pbCntr.uPncl,pbCntr.dvPncl,stats);
				

			}
                        
                        
                        //Compute (u.grad)u on the bundle.
			pbCntr.dvPncl.fillPencil(0.0);
                        udotgradu(pbCntr.uBndl,pbCntr.dvPncl,xfac,yfac,zfac);
			//pbCntr.rhskPncl = pbCntr.rhskPncl + pbCntr.dvPncl;
			pbCntr.rhskPncl += pbCntr.dvPncl;
			
                        //Compute vector laplacian on the bundle.
			pbCntr.dvPncl.fillPencil(0.0);
			vlapl(pbCntr.uBndl,pbCntr.dvPncl,xfac*xfac,yfac*yfac,zfac*zfac);
                        //pbCntr.dvPncl = pbCntr.dvPncl*(1.0/rho);
			//pbCntr.rhskPncl = pbCntr.rhskPncl + pbCntr.dvPncl;
			pbCntr.rhskPncl += pbCntr.dvPncl*(1.0/rho);
			
                        //Compute force pencil
			pbCntr.fPncl.fillPencil(0.0);
                        set_force(pbCntr.fPncl,params,grid,j);
			//pbCntr.rhskPncl = pbCntr.rhskPncl + pbCntr.fPncl;
			pbCntr.rhskPncl += pbCntr.fPncl;
			
                        //Copy from pencil to mesh
                        pencil2ff(pbCntr.rhskPncl,meshCntr.RHSk,j,k);
			
                }
        }
	//std::cout << "umax after calc: " << stats.umax << std::endl;
	
        
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
void calc_gradpsi(MeshContainer &meshCntr, PBContainer &pbCntr, const SolverParams &params, const Grid &grid)
{
        size_t Ny = meshCntr.psi.ny_;
        size_t Nz = meshCntr.psi.nz_;

        Real xfac=grid.invdx,yfac=grid.invdy,zfac=grid.invdz;

        // Loop over mesh (y,z)-plane 
        for (size_t j=0;j<Ny;j++)
        {
                for (size_t k=0;k<Nz;k++)
                {
                        ff2bundle(meshCntr.psi,pbCntr.psiBndl,j,k); //Copy
                        sgrad(pbCntr.psiBndl,pbCntr.dvPncl,xfac,yfac,zfac); //Compute gradient on bundle
                        pencil2ff(pbCntr.dvPncl,meshCntr.gradpsi,j,k);
                        
                }
        }
}

/// Enforces solenoidal condition on velocity field by subtracting grad(psi) from ustar.
void enforce_solenoidal(MeshContainer &meshCntr, PBContainer &pbCntr, const SolverParams &params, const Grid &grid, const Int k_rk)
{
        ///Set param constants
        Real dt = params.dt, rho = params.rho;
        
        /// Create RK3 coefficient
        RK3Coeff coeffs;
        Real alphak = coeffs.alpha(k_rk);

        apply_pbc(meshCntr.psi);
        calc_gradpsi(meshCntr,pbCntr,params,grid);
        meshCntr.u = meshCntr.ustar-(meshCntr.gradpsi*(2*alphak*dt/rho));
}

void update_timestep(SolverParams &params, Stats &stats, const Grid &grid)
{
	if (params.currentTimestep !=1) {
		Real c1=0.1;
		Real c2=c1; //Courant numbers for advection and diffusion respectively

		Real nu = params.viscosity,dx = grid.dx,L=grid.dx;
        
		//Factor of 1/3 since dx=dy=dz
	
		Real adv = (1.0/3)*c1*dx/stats.umax;
		Real diff = (1.0/3)*c2*pow(L,2)/nu;
		//std::cout << "adv : " << adv << "\t diff: " << diff << std::endl;
		//Set new time step size according to CFL condition
		params.dt = std::min(adv,diff);
		std::cout << params.dt << std::endl;
	}

}

/// Function that performs RK3 substeps from k=1 to k=3
void RK3_stepping(MeshContainer &meshCntr, PBContainer &pbCntr, SolverParams &params, Stats &stats, const Grid &grid)
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
                calc_RHSk(meshCntr,pbCntr,params,stats,grid,k_rk); //diagnostics also calculated here
                // If k_rk == 1 update the timestep dt
                if (k_rk == 1)
                        update_timestep(params,stats,grid);
                
                calc_ustar(meshCntr,pbCntr,params,grid,k_rk);

                /// Solve the Poisson eq for Psi.
                apply_pbc(meshCntr.ustar);
                calc_divustar(meshCntr,pbCntr,params,grid,k_rk);
                PoissonEq(meshCntr,params,grid,k_rk);

                /// Update pressure and enforce solenoidal condition
                update_pressure(meshCntr);
                enforce_solenoidal(meshCntr,pbCntr,params,grid,k_rk);
        }
        
} //End RK3_stepping()

void save_data_at_step(const MeshContainer &meshCntr, const SolverParams &params, const Stats &stats, const Grid &grid, const Int t)
{
	std::string kfname = kolmofname("kolmo",params.kf);

	std::string statsfname = step_fname(kfname, ".stats",meshCntr.u.nx_,t);

	for (Int i=0;i<3;i++)
	{
		std::string meshfname = step_fname(kfname+"_component_"+std::to_string(i), ".dat",meshCntr.u.nx_,t);
		writeToFile_1DArr(meshfname,meshCntr.u,i,grid); //Write i-component of u field to file
	}
	
	writeStatsToFile(statsfname,meshCntr,stats,grid,t); //Calculate and write statistics to file
	//std::cout << "wrote " << meshfname << " and " << statsfname << " to file at timestep " << t << std::endl;
}

void execprint(const MeshContainer &meshCntr, const SolverParams &params, const Stats &stats, const Grid &grid)
{
	std::cout << "Started run with N: " << meshCntr.u.nx_ << " maxsteps: " << params.maxTimesteps <<
		", saving every " << params.saveintrvl << " timesteps." << std::endl;
	std::cout << std::endl;
}

Int main()
{
        auto t1 = Clock::now();
 
        /// Time settings.
        const Int maxtsteps = 10;

        /// Create parameter object and initialise parameters.
        SolverParams params;
        params.maxTimesteps = maxtsteps;
        params.currentTimestep = 0;
        params.kf = 1.0; //Kolmogorov frequency
        params.rho = 1.0;
        params.viscosity = 1.0;
	params.saveintrvl = 4;
	
        
        /// Set grid sizes
        const Real L0 = -M_PI, L1 = M_PI; // x,y,z in [-pi,pi]
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

	/// Create Bundle & Pencil objects
	Bundle uBndl(Nsize,3);
	Bundle ustarBndl(Nsize,3);
	Bundle pBndl(Nsize,1);
	Bundle psiBndl(Nsize,1);
	Pencil uPncl(Nsize,3);
	Pencil ustarPncl(Nsize,3);
	Pencil dvPncl(Nsize,3);
	Pencil dsPncl(Nsize,1);
	Pencil forcePncl(Nsize,3);
	Pencil rhskPncl(Nsize,3);
	Pencil rhsk_1Pncl(Nsize,3);
	
        /// Create container object with references to meshes
        MeshContainer meshCntr(uu,ustar,RHSk,RHSk_1,pp,psi,gradpsi,fftwPsi);
	PBContainer pbCntr(uBndl,ustarBndl,pBndl,psiBndl,uPncl,ustarPncl,dvPncl,dsPncl,
			   forcePncl,rhskPncl,rhsk_1Pncl);
	
        /// At first we get initial timestep size from forcing
        Pencil initforcePncl(Nsize,1);
        for (size_t i=0;i<uu.nx_;i++)
                initforcePncl(i,0,0) = sin(params.kf*grid.x(i));

	/// Initial set of the timesteps
        stats.calc_pncl_absmax(initforcePncl);
	//std::cout << stats.umax << std::endl;
        stats.umax_old = stats.umax;
        update_timestep(params,stats,grid);
        stats.umax = 0;
        //std::cout << 0 << " : " << params.dt << std::endl;
        //Set umax to zero again to get actual calculation from run
	
/// Run Runge-Kutta 3 substepping
	execprint(meshCntr,params,stats,grid);
	
        for (Int ts = 1;ts<params.maxTimesteps;ts++)
        {
		meshCntr.u.print();
		params.currentTimestep = ts;
		
                RK3_stepping(meshCntr,pbCntr,params,stats,grid);
		
		//Save data every 4th timestep
		//if (ts % params.saveintrvl == 0)
			//save_data_at_step(meshCntr,params,stats,grid,ts);
			
		//calc_mesh_avg(meshCntr,stats,grid);
                //stats.urms = stats.calc_mesh_rms(meshCntr.u);
                //calc_Re(params,stats,grid);
                //std::cout << ts << " : " << params.dt << std::endl;
	}
	meshCntr.u.print();
	
	//save last step again to be sure
	//save_data_at_step(meshCntr,params,stats,grid,params.maxTimesteps+1);
	//apply_pbc(meshCntr.u);
	//calc_curlu(meshCntr,curlu,pbCntr,params,stats,grid);
        //writeCurl(curlu,params,stats,grid,params.maxTimesteps+1);
        auto t2 = Clock::now();
        std::cout << "Delta t2-t1: " 
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds" << std::endl;
        
        return 0;
        
} //End main()
