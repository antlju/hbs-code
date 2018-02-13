#include "includes.h"

/// Test program for computing the curl of the ABC-field.
/// - version 1.

/// Takes two vector fields and produces a scalar function.
void dotprod(const Mesh &A, const Mesh &B, Mesh &C)
{

	for (size_t i=0;i<A.nx_;i++)
	{
		for (size_t j=0;j<A.ny_;j++)
		{
			for (size_t k=0;k<A.nz_;k++)
			{
                                
				C(i,j,k) = A(i,j,k,0)*B(i,j,k,0)
					+A(i,j,k,1)*B(i,j,k,1)
					+A(i,j,k,2)*B(i,j,k,2);
			}
		}
	}

}

void setup(Mesh &u, const Pencil &x, const Pencil &y, const Pencil &z)
{
	//Real A=1.0,B=1.0,C=1.0;
	Real L=1.0;
	for (size_t i=0;i<u.nx_;i++)
	{
		for (size_t j=0;j<u.ny_;j++)
		{
			for (size_t k=0;k<u.nz_;k++)
			{
				u(i,j,k,0) = exp(L*x(i))+exp(L*y(j))+exp(L*z(k));
				//u(i,j,k,0) = pow(sin(x(i)),2)+pow(sin(y(j)),2)+pow(sin(z(k)),2);
				//u(i,j,k,1) = sin(y(j));
				//u(i,j,k,2) = sin(z(k));
				//u(i,j,k,0) = x[i];
				//u(i,j,k,1) = y[i];
				//u(i,j,k,2) = z[i];
			}
		}
	}

}

void setup2(Mesh &u, const Pencil &x, const Pencil &y, const Pencil &z)
{
	//Real A=1.0,B=1.0,C=1.0;
	Real L=1.0;
	for (size_t i=0;i<u.nx_;i++)
	{
		for (size_t j=0;j<u.ny_;j++)
		{
			for (size_t k=0;k<u.nz_;k++)
			{
				u(i,j,k,0) = exp(L*x(i))+exp(L*y(j))+exp(L*z(k));
				//u(i,j,k,0) = -2*(sin(x(i))+sin(y(j))+sin(z(k)));
				//u(i,j,k,0) = x[i];
				//u(i,j,k,1) = y[i];
				//u(i,j,k,2) = z[i];
			}
		}
	}

}

void apply_lapl(const Mesh &u, Mesh &divu, const Real xfac, const Real yfac, const Real zfac)
{
        Bundle Dbundle(u.nx_,1); /// scalar bundle
        Pencil Dpencil(u.nx_,1); /// scalar pencil
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {
                        ff2bundle(u,Dbundle,j,k);
                        lapl(Dbundle,Dpencil,xfac,yfac,zfac);
                        pencil2ff(Dpencil,divu,j,k);
                }
        }
}

Int main()
{
        Int Nx = NX;
        Int Ny = NY; 
        Int Nz = NZ;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        Pencil x(Nx);
        x.printpencil();
        linspace(x,L0,L1,dx);
        x.printpencil();

        Mesh u(Nx,Ny,Nz,1);
        setup(u,x,x,x);
        //apply_pbc(u);

	Mesh mangu(Nx,Ny,Nz,1);
		
	Mesh g(Nx,Ny,Nz,1);
	setup2(mangu,x,x,x);
	apply_lapl(u,g,1.0/(dx*dx),1.0/(dx*dx),1.0/(dx*dx));

	//printfield(u);
	//printfield(g);
	
	for (size_t vii=0;vii<g.nvar_;vii++)
	{
		std::cout << "------- component: " << vii << " -------" << std::endl;
		std::cout << "---------------------------------------" << std::endl;
		for (size_t ii=0;ii<g.nx_;ii++)
		{
			for (size_t jj=0;jj<g.ny_;jj++)
			{
				for (size_t kk=0;kk<g.nz_;kk++)
				{
					//std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii,jj,kk,vii) << std::endl;
					if (mangu(ii,jj,kk,vii) < 1e-15)
					{
						std::cout << "div by 0" << "\t";
					}
					else
					{
						
						std::cout << g(ii,jj,kk,vii)/u(ii,jj,kk,vii) << "\t";
							
					}
						

				}
				std::cout << "---------------------------------------" << std::endl;
				std::cout << std::endl;
			}
		}
	}
        
        
        
	return 0;
}

