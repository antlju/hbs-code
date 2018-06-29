#include "includes.h"
#include <chrono>

using Clock=std::chrono::high_resolution_clock;

void init_sinz(Mesh &u, Grid &grid)
{
	Int NN = u.nx_;
	Pencil z = grid.x;
	for (Int i=0;i<NN;i++)
	{
		for (Int j=0;j<NN;j++)
		{
			for (Int k=0;k<NN;k++)
			{
				//h[fIdx(i,j,k)] = i+NX*(j+NY*k);
				u(i,j,k,0) = sin(z(k));
			}
		}
	}
}

void calc_zderiv(const Mesh& u, Mesh &du, const Grid &grid)
{
	size_t Ny = u.ny_;
	size_t Nz = u.nz_;
	size_t Nx = u.nx_;
	Real zfac = grid.invdz;

	Bundle B(Nx,1);
	Pencil P(Nx,1);
	for (size_t j=0;j<Ny;j++)
	{
		for (size_t k=0;k<Nz;k++)
		{
			ff2bundle(u,B,j,k);
			for (size_t i=0;i<Nx;i++)
			{
				P(i) = delz(B,zfac,i,0);
			}

			pencil2ff(P,du,j,k);
			
		}
	}
}

Real calcMaxErr(const Mesh &du, const Grid &grid)
{
	Real maxE = 0;
	for (size_t i=0;i<du.nx_;i++)
	{
		for (size_t j=0;j<du.ny_;j++)
		{
			for (size_t k=0;k<du.nz_;k++)
			{
				Real val = fabs(du(i,j,k)-cos(grid.x(k)));
				if (val > maxE)
					maxE = val;
			}
		}
	}
	return maxE;
}

Int main()
{
	
	const Int NN = 256;
	/// Set grid sizes
        const Real L0 = 0, L1 = 2*M_PI; // x,y,z in [0,2pi]

	Grid grid(NN,NN,NN,L0,L1);

	Mesh uu(NN,NN,NN,1);
	Mesh duu(NN,NN,NN,1);

	init_sinz(uu,grid);
	apply_pbc(uu);

	auto t1 = Clock::now();
	calc_zderiv(uu,duu,grid);
	auto t2 = Clock::now();

	//duu.print();
	Real maxErr = calcMaxErr(duu,grid);
	std::cout << "Max error: " << maxErr << std::endl;
	std::cout << "Time taken for CPU z-deriv: "
		  << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
		  << " ms" << std::endl;
	std::cout << "System size (N = " << NN <<")^3" << std::endl;
	
}
