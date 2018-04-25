#include "includes.h"

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
                                Real term = sin(x(i))+sin(x(j))+sin(x(k));
                                u(i,j,k,0) = term;
                                u(i,j,k,1) = term;
                                u(i,j,k,2) = term;
                        }
                }
        }
               

}


/// Function for copying from a bundle to a pencil.
void bundle2pencil(const Bundle &B, Pencil &P)
{
        assert((B.nx_ == P.nx_) && (B.nvars_ == P.nvars_) );
        for (size_t i=0;i<B.nx_;i++)
        {
                for (size_t vi=0;vi<B.nvars_;vi++)
                {
                        P(i,0,vi) = B(i,0,vi);
                }
        }
}


void set_force(Pencil &force, const Pencil &y, const Int j)
{
        for (size_t i=0;i<force.nx_;i++)
        {
                force(i,0,0) = sin(y(j));
                force(i,0,1) = 0.0;
                force(i,0,2) = 0.0;
        }
}

void compute_ustar(const Mesh &ff, Mesh &dff, Mesh &pp,const Real nu, const Real xfac, const Real yfac, const Real zfac, const Pencil &y)
{
        Bundle u(ff.nx_,3); /// Vector bundle
        Pencil ustar(ff.nx_,3); /// Vector pencil
        Pencil uDelu(ff.nx_,3); /// Vector pencil
        Pencil Laplu(ff.nx_,3); /// Vector pencil
        Pencil force(ff.nx_,3); /// Vector pencil
        Pencil uPencil(ff.nx_,3);
        
        for (size_t j=0;j<ff.ny_;j++)
        {
                for (size_t k=0;k<ff.nz_;k++)
                {
                        ff2bundle(ff,u,j,k); //Copy
                        
                        bundle2pencil(u,uPencil); //We need the old u in computing u*.
                        
                        udotgradu(u,uDelu,xfac,yfac,zfac); //Compute (u.grad)u on the bundle.
                        
                        vlapl(u,Laplu,xfac*xfac,yfac*yfac,zfac*zfac); //Compute vector laplacian on the bundle.
                        
                        set_force(force,y,j); //Compute the correct forcing term from Kolmogorv flow along the current (j,k) pencil.
                        
                        ustar = uPencil + uDelu + Laplu*nu + force; //Set u*-pencil
                        //ustar = up + uDelu + Laplu + force;
                        
                        pencil2ff(ustar,dff,j,k);
                        
                }
        }

        
        
}

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

        
        
}

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

}
Int main()
{
        /// Init dimensions.
        Int Nx=NX,Ny=NY,Nz=NZ;

        /// Define computation grid
        Real L0=-2*M_PI;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx-1);
        Real xlen = (L1-L0);
        Pencil x(Nx);
        linspace(x,L0,L1,dx);

        /// Create objects for calculation and allocate appropriate space as.
        Mesh ff(Nx,Ny,Nz,3);   //main velocity vector field
        Mesh dff(Nx,Ny,Nz,3); //space for computations
        Mesh pp(Nx,Ny,Nz,1); //pressure scalar field
        Mesh psi(Nx,Ny,Nz,1); //result from Poisson equation.
        Mesh gradpsi(Nx,Ny,Nz,3); //Gradient of psi to enforce solenoidal condition.
        fftwMesh psifftw(Nx,Ny,Nz); //Space for 3D FFTW solver for Poisson eq.

        /// Set initial condition and apply periodic boundaries.
        Real Re = sqrt(2.0)-0.00001; // Reynolds number
        init_ff(ff,x);
        apply_pbc(ff);

        /// Compute u*
        compute_ustar(ff,dff,pp,1.0/Re,1.0/dx,1.0/dx,1.0/dx,x);

        /// Compute div(u*), the source term for the Poisson equation.
        compute_divustar(dff,psi,1.0/dx,1.0/dx,1.0/dx);

        /// Solve the Poisson equation Lapl(psi) = div(ustar)
        PoissonSolve3D(psi,psifftw,0,xlen);
        
        // Divide by FFTW backtransform scaling.
        psi = psi*(1.0/(Nx*Ny*Nz));

        // Update pressure
        pp = pp + psi;

        //Compute gradient of psi
        compute_gradpsi(psi,gradpsi,1.0/dx,1.0/dx,1.0/dx);
        //Enforce solenoidal condition on u
        dff = dff-gradpsi;
}
