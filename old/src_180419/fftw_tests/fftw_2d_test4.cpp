#include "includes.h"

#include <fftw3.h>

void initialise(Real *in, Real q1, Real q2,Pencil &x,Pencil &y)
{
	size_t Nx = x.nx_,Ny = y.nx_;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        Real val = sin(q1*x(i))*sin(q2*y(j));
                        in[j+(Nx+2)*i] = cos(q1*x(i))*cos(q2*y(j));
		}
	}
}
	       
Int main()
{
        Real q1=1.0,q2=1.0;
        Real pi = M_PI;
        Int Size = 16;
	Int Nx=Size,Ny=Size; 
	fftw_complex *mem_;
	mem_ = (fftw_complex*)
			fftw_malloc(
                                sizeof(fftw_complex) * (Nx*(2*(Ny/2+1)))
				);
	
	//fftw_complex *out = mem_;
	Real *in = mem_[0];
	fftw_complex *out = mem_;
	fftw_plan forward = fftw_plan_dft_r2c_2d(Nx,Ny,
                                                 in,out,
                                                 FFTW_MEASURE);
        fftw_plan backward = fftw_plan_dft_c2r_2d(Nx,Ny,
                                                  out,in,
                                                  FFTW_MEASURE);

	Real L0=0;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx);
	Pencil x(Nx);
        linspace(x,L0,L1,dx);
        Real dxx = x(1)-x(0);
        
	
	initialise(in,q1,q2,x,x);
        

	fftw_execute(forward);

        for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<(Ny/2)+1;j++)
		{
                        Real InvSize = L1-L0;
                        Int JJ = 0,II=i;
                        //if (2*i<Nx)
                        //        II = i;
                        //else
                        //        II = Nx-i;
                        
                        if (2*i<Nx)
                                II = i;
                        else
                                II = Nx-i;

                        Real fac = (pow(2*pi*II/InvSize,2)+pow(2*pi*j/InvSize,2));
                        if (fac == 0)
                        {
                                out[j+Nx*i][0] = 0;
                                out[j+Nx*i][1] = 0;
                        }
                        else
                        {
                                out[j+Nx*i][0] = (-1.0/fac)*out[j+Nx*i][0];
                                out[j+Nx*i][1] = (-1.0/fac)*out[j+Nx*i][1];
                        }
                }
        }

        
        fftw_execute(backward);

        //std::cout.precision(4);
        //Printing result of inverse transform compared to analytical
        for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        //Real ana_psi = -sin(q1*x(i))*sin(q2*x(j))/(pow(q1,2)+pow(q2,2));
                        //Real ana_A = sin(q1*x(i))*sin(q2*x(j));
                        //Real diff = fabs(in[j+(Nx+2)*i]/(Nx*Ny)-ana_psi);
                        //if (diff < 1e-14)
                        //        diff = 0.0;

                        
                        std::cout << x(i) <<"\t"<< x(j) << "\t" << in[j+(Nx+2)*i]/(Nx*Ny) << std::endl;
                        
		}
	}
        
	fftw_destroy_plan(forward);
        fftw_destroy_plan(backward);
	fftw_free(mem_);








	
        return 0;
}
