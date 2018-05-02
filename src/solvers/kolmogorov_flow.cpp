#include "includes.h"
std::string set_fname(std::string base, std::string end, Int step)
{
        std::string timeval = "_t_";
        std::string timeNo = "_tNum_";
        return base+timeNo+std::to_string(step)+end;
}

Int writeToFile(std::string fname, const Mesh &u, const Int vi, const Pencil &x, const Pencil &y, const Pencil &z)
{
        std::ofstream openfile("/home/anton/dev/hbs/simdata/"+fname, std::ios::trunc);
        openfile << "i,j,k,x,y,z,f" << std::endl;
        //for (size_t i=0;i<u.nx_;i++)
        //{
        openfile << "#---------------------------------" << std::endl;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                /*
                                  openfile << i << "," << j << "," << k << ","
                                  << x(i) << "," << y(j) << "," << z(k) << ","
                                  << u(i,j,k,vi) << std::endl;
                                */

                                openfile << u(i,j,k,vi) << "\t";
                        }
                        openfile << std::endl;       
                }
                openfile << "#---------------------------------" << std::endl;
        }
                // }
        openfile.close();
        std::cout << "Wrote " << fname << std::endl;
        return 0;

}

/// Set initial conditions.
void init_ff(Mesh &u, const Pencil &x)
{
        //size_t Ng = u.ng_,Nvars=u.nvar_;
        size_t Nx=u.nx_,Ny=u.ny_,Nz=u.nz_;
        
        for (size_t i=0;i<Nx;i++)
        {
                for (size_t j=0;j<Ny;j++)
                {
                        for (size_t k=0;k<Nz;k++)
                        {
                                //Real term = sin(x(i))+sin(x(j))+sin(x(k));
                                Real term = sin(x(j));
                                u(i,j,k,0) = 0.0;
                                u(i,j,k,1) = 0.0;
                                u(i,j,k,2) = 0.0;
                        }
                }
        }
               

}


void set_force(Pencil &force, const Pencil &y, const Int j)
{
        Real Amp = 1.0;
        for (size_t i=0;i<force.nx_;i++)
        {
                force(i,0,0) = Amp*sin(y(j));
                force(i,0,1) = 0.0;
                force(i,0,2) = 0.0;
        }
        
} //End set_force

void compute_ustar(const Mesh &ff, Mesh &dff, Mesh &pp,
                   const Real nu, const Real dt, const Real xfac, const Real yfac, const Real zfac,
                   const Pencil &y)
{
        
        Bundle u(ff.nx_,3); /// Vector bundle
        Bundle pBundle(pp.nx_,1); /// Scalar bundle
        Pencil ustar(ff.nx_,3); /// Vector pencil
        Pencil uDelu(ff.nx_,3); /// Vector pencil
        Pencil Laplu(ff.nx_,3); /// Vector pencil
        Pencil force(ff.nx_,3); /// Vector pencil
        Pencil uPencil(ff.nx_,3);
        Pencil pPencil(ff.nx_,3); // Vector Pencil
        for (size_t j=0;j<ff.ny_;j++)
        {
                for (size_t k=0;k<ff.nz_;k++)
                {
                        ff2bundle(ff,u,j,k); //Copy velocity and pressure
                        ff2bundle(pp,pBundle,j,k);
                        
                        bundle2pencil(u,uPencil); //We need the old u in computing u*.
                        
                        udotgradu(u,uDelu,xfac,yfac,zfac); //Compute (u.grad)u on the bundle.
                        
                        vlapl(u,Laplu,xfac*xfac,yfac*yfac,zfac*zfac); //Compute vector laplacian on the bundle.

                        sgrad(pBundle,pPencil,xfac,yfac,zfac);
                        
                        set_force(force,y,j); //Compute the correct forcing term from Kolmogorv flow along the current (j,k) pencil.
                        
                        ustar = uPencil - pPencil*dt - uDelu*dt + Laplu*(nu*dt) + force*dt; //Set u*-pencil
                        //ustar = up + uDelu + Laplu + force;
                        
                        pencil2ff(ustar,dff,j,k);
                        
                }
        }

        
        
} //End compute_ustar()

void compute_divustar(const Mesh &ustar, Mesh &divustar, const Real xfac, const Real yfac, const Real zfac)
{
        Bundle ustarB(ustar.nx_,3); /// Vector bundle
        Pencil divustarP(divustar.nx_,1); /// Scalar pencil

        for (size_t j=0;j<ustar.ny_;j++)
        {
                for (size_t k=0;k<ustar.nz_;k++)
                {
                        ff2bundle(ustar,ustarB,j,k); //Copy
                        div(ustarB,divustarP,xfac,yfac,zfac); // Compute divergence on bundle
                        pencil2ff(divustarP,divustar,j,k);
                        
                }
        }

        
        
} //End compute_divustar()

void compute_gradpsi(const Mesh &psi, Mesh &gradpsi,const Real xfac, const Real yfac, const Real zfac)
{
        Bundle B(psi.nx_,psi.nvar_); //Scalar bundle
        Pencil P(psi.nx_,3); // Vector pencil

        for (size_t j=0;j<psi.ny_;j++)
        {
                for (size_t k=0;k<psi.nz_;k++)
                {
                        ff2bundle(psi,B,j,k); //Copy
                        sgrad(B,P,xfac,yfac,zfac); //Compute gradient on bundle
                        pencil2ff(P,gradpsi,j,k);
                        
                }
        }

} //End compute_gradpsi()

void compute_rhs(const Mesh &ff, Mesh &dff, Mesh &pp, Mesh &psi, Mesh &gradpsi, fftwMesh &psifftw,
                 const Real dt, const Real dx, const Real dy, const Real dz,
                 const Real xlen, const Real Re, const Pencil &y)
{
        /// Compute u*
        compute_ustar(ff,dff,pp,1.0/Re,dt,1.0/dx,1.0/dy,1.0/dz,y);
        apply_pbc(dff);
        
        /// Compute div(u*), the source term for the Poisson equation.
        /// We save div(u*) in psi to save space of one scalar field.
        compute_divustar(dff,psi,1.0/dx,1.0/dy,1.0/dz);

        /// Solve the Poisson equation Lapl(psi) = div(ustar)
        PoissonSolve3D(psi,psifftw,0,xlen);
        
        // Divide by FFTW backtransform scaling.
        psi = psi*(1.0/(ff.nx_*ff.ny_*ff.nz_));
        apply_pbc(psi);
        
        // Update pressure
        pp = pp + psi;

        //Compute gradient of psi
        compute_gradpsi(psi,gradpsi,1.0/dx,1.0/dx,1.0/dx);
        
        //Enforce solenoidal condition on u
        dff = dff-(gradpsi*(2*dt));
        apply_pbc(dff);
        
} //End compute_rhs()

void RK4(Mesh &ff, Mesh &dff, Mesh &pp, Mesh &psi, Mesh &gradpsi, fftwMesh &psifftw,
         Mesh &ff_step, Mesh &k1, Mesh &k2, Mesh &k3, Mesh &k4,
         const Real dt, const Int maxtsteps,
         const Real dx, const Real dy, const Real dz, const Real xlen, const Real Re,
         const Pencil &y)
{
        Real c1=1.0/6.0,c2=1.0/3.0;

        for (Int tstep=1;tstep<maxtsteps;tstep++)
        {
                compute_rhs(ff,dff,pp,psi,gradpsi,psifftw,
                            dt,dx,dy,dz,
                            xlen,Re,y);
                /*
                k1 = dff*dt;

                ff_step = ff+(k1*0.5);
                apply_pbc(ff_step);
                compute_rhs(ff_step,dff,pp,psi,gradpsi,psifftw,
                            dt,dx,dy,dz,
                            xlen,Re,y);
                k2 = dff*dt;

                ff_step = ff+(k2*0.5);
                apply_pbc(ff_step);
                compute_rhs(ff_step,dff,pp,psi,gradpsi,psifftw,dt,
                            dx,dy,dz,
                            xlen,Re,y);
                k3 = dff*dt;

                ff_step = ff+(k3*0.5);
                apply_pbc(ff_step);
                compute_rhs(ff_step,dff,pp,psi,gradpsi,psifftw,
                            dt,dx,dy,dz,
                            xlen,Re,y);
                k4 = dff*dt;

                ff = ff + k1*c1 + k2*c2 + k3*c2 + k4*c1;
                */
                ff = ff+dff*dt;
                apply_pbc(ff);
                // std::cout << "Completed RK4 step #"<< tstep << std::endl;
                writeToFile(set_fname("kolmo_v1_component_0_",".csv",tstep),ff,0,y,y,y);
                writeToFile(set_fname("kolmo_v1_component_1_",".csv",tstep),ff,1,y,y,y);
                writeToFile(set_fname("kolmo_v1_component_2_",".csv",tstep),ff,2,y,y,y);
        }
        
} //End RK4()


Int main()
{
        /// Time settings.
        Real dt = 0.01;
        Int maxtsteps = 100;
        
        /// Init dimensions.
        Int NN= 4;
        Int Nx=NN,Ny=NN,Nz=NN;

        /// Define computation grid
        Real L0=-2*M_PI;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx-1);
        Real xlen = (L1-L0);
        Pencil x(Nx);
        linspace(x,L0,L1,dx);

        /// Create objects for calculation and allocate appropriate space as.
        Mesh ff(Nx,Ny,Nz,3);   //main velocity vector field
        Mesh dff(Nx,Ny,Nz,3); //space for computations on velocity field
        Mesh pp(Nx,Ny,Nz,1); //pressure scalar field
        Mesh psi(Nx,Ny,Nz,1); //result from Poisson equation.
        Mesh gradpsi(Nx,Ny,Nz,3); //Gradient of psi to enforce solenoidal condition.
        fftwMesh psifftw(Nx,Ny,Nz); //Space for 3D FFTW solver for Poisson eq.

        /// Runge Kutta 4 step objects
        Mesh k1(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k2(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k3(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k4(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh ff_step(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        
        /// Set initial condition and apply periodic boundaries.
        Real Re = sqrt(2.0)-0.1; // Reynolds numbre
        init_ff(ff,x);
        apply_pbc(ff);
        ff.print();
        
        /// Run solver with RK4 timestepping.
        RK4(ff,dff,pp,psi,gradpsi,psifftw,
            ff_step,k1,k2,k3,k4,
            dt,maxtsteps,
            dx,dx,dx,xlen,Re,
            x);

        ff.print();
}
