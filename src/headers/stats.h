#pragma once

#include "vartypedef.h"
#include "typedef.h"
#include <cmath>

/// Statistics calculation and storage class

class Stats {
public:
        /// max(abs(u)), The maximum absolute velocity
        Real umax;
        Real umax_old;

        Real urms;

        /// Function for calculating statistics on pencils
        /// ---------------------------------------------------


        
        /// ---------------------------------------------------

        
        /// Functions for calculating statistics on mesh
        /// -----------------------------------------------------
        Real calc_mesh_max(const Mesh &inmesh)
        {
                Real maxval = inmesh(0,0,0);
                for (size_t i=0;i<inmesh.nx_;i++)
                {
                        for (size_t j=0;j<inmesh.ny_;j++)
                        {
                                for (size_t k=0;k<inmesh.nz_;k++)
                                {
                                        if (inmesh(i,j,k) >= maxval)
                                                maxval = (Real)inmesh(i,j,k);
                                }
                        }
                }

                return maxval;
        }

        Real calc_mesh_rms(const Mesh &inmesh)
        {
                Real sum;
                for (size_t i=0;i<inmesh.nx_;i++)
                {
                        for (size_t j=0;j<inmesh.ny_;j++)
                        {
                                for (size_t k=0;k<inmesh.nz_;k++)
                                {
                                        sum += pow(inmesh(i,j,k),2);
                                }
                        }
                }
                Int Nx = inmesh.nx_,Ny=inmesh.ny_,Nz=inmesh.nz_;
                Real ret = sqrt(sum/(Nx*Ny*Nz));
                
                return ret;
        }

        /// ---------------------------------------------------------
            
}; //End class Stats
