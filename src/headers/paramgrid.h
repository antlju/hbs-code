#pragma once

#include <vector>
#include <iostream>
#include "vartypedef.h"
#include "pmesh.h"
#include "typedef.h"

/// A parameter container class for use in various solvers.
class SolverParams {
public:
        /// Time step size (to be updated according to CFL before each Runge-Kutta ingegration step).
        Real dt; 
        Real dt_old; /// Previous step size
        Real dt_init;
        
        /// Max time steps.
        Int maxTimesteps;

        ///current time step
        Int currentTimestep;
        
        /// Reynolds number, viscosity
        /// Might require updates in time steps.
        Real Re;
        Real viscosity;

        /// Density
        Real rho;

        /// Kolmogorov frequency
        Real kf;

}; // End class SolverParams


/// A class for finite difference grid properties
class Grid {
public:
        /// Grid sizes
        size_t nx_,ny_,nz_;
        
        /// Spatial step sizes
        Real dx,dy,dz;
        Real invdx,invdy,invdz;
        Real invdx2,invdy2,invdz2;

        /// Start- and endpoints
        Real L0_,L1_; // For uniform grids
        Real L0x,L1x;
        Real L0y,L1y;
        Real L0z,L1z;

        /// Axes (most often as linspaces)
        Pencil x;
        Pencil y;
        Pencil z;

        /// Total length of axes
        Real xlen,ylen,zlen;

        void xlinspace()
        {
                for (size_t i=0;i<x.nx_;i++)
                {
                        x(i) = L0x+dx*i;
                        //std::cout << u(i) << std::endl;
                }
        }

        /// Constructor for uniform grid, i. e same size along all three dimensions.
        Grid(Int Nx,Int Ny,Int Nz,Real L0, Real L1) :
                /// Init dims
                nx_(Nx),ny_(Ny),nz_(Nz),
                /// Init start and endpts
                L0_(L0),L1_(L1),L0x(L0),L1x(L1),L0y(L0),L1y(L1),L0z(L0),L1z(L1),
                /// Init pencils
                x(Nx),y(Ny),z(Nz)
        {
                /// Calculate step size and lengths
                dx = (L1-L0)/(Nx-1); dy = dx; dz = dx;
                xlen = (L1-L0); ylen = xlen; zlen = xlen;

                /// Create x-axis linspace.
                xlinspace();

                /// Inverse steps for use in differential operators
                invdx = 1.0/dx; invdy = 1.0/dy; invdz = 1.0/dz;
                invdx2 = pow(invdx,2); invdy2 = pow(invdy,2); invdz2 = pow(invdz,2);
        }

                
}; // End class Grid

/// Container class for Runge-Kutta 3 coefficients 
class RK3Coeff {
public:

        /// Coefficient values taken from Rosti & Brandt 2017
        Real alpha1 = 4.0/15;
        Real alpha2 = 1.0/15;
        Real alpha3 = 1.0/6;
        
        Real beta1 = 8.0/15;
        Real beta2 = 5.0/12;
        Real beta3 = 3.0/4;
        
        Real gamma1 = 0.0;
        Real gamma2 = -17.0/60;
        Real gamma3 = -5.0/12;

        /// Functions for returning constants for each k
        Real alpha(Int k)
        {
                assert( k >= 1 && k <= 3);
                
                switch(k)
                {
                case 1: return alpha1;
                case 2: return alpha2;
                case 3: return alpha3;
                default: return 0.0;
                }
                
        }
        
        Real beta(Int k)
        {
                assert( k >= 1 && k <= 3);
                
                switch(k)
                {
                case 1: return beta1;
                case 2: return beta2;
                case 3: return beta3;
                default: return 0.0;
                }
                
        }
        
        Real gamma(Int k)
        {
                assert( k >= 1 && k <= 3);
                
                switch(k)
                {
                case 1: return gamma1;
                case 2: return gamma2;
                case 3: return gamma3;
                default: return 0.0;
                }
                
        }
        
};
