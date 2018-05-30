
#pragma once

#include "vartypedef.h"
#include "typedef.h"
#include <cmath>

/// Statistics calculation and storage class

class Stats {
public:
        /// max(abs(u)), The maximum absolute velocity
	Real Re;
	
        Real umax;
	
	Real energy;
	
	Real urms;

	Int isSaveStep = 0;

	        
        /// ---------------------------------------------------
        /// Functions for calculating statistics on mesh
	/// -----------------------------------------------------
	
	void calc_mesh_urms(const Mesh &u)
	{
		urms = 0.0;
		//for (size_t vi=0;vi<u.nvar_;vi++)
		for (size_t vi=0;vi<1;vi++) // only x-component
		{
			for (size_t i=0;i<u.nx_;i++)
			{
				for (size_t j=0;j<u.ny_;j++)
				{
					for (size_t k=0;k<u.nz_;k++)
					{
						urms += u(i,j,k,vi)*u(i,j,k,vi);
					}
				}
			}
		}
		//urms = sqrt(urms/(u.nvar_*u.nx_*u.ny_*u.nz_));
		urms = sqrt(urms/(u.nx_*u.ny_*u.nz_));
	}

	void calc_mesh_energy(const Mesh &inmesh)
	{
		//for (size_t vi=0;vi<inmesh.nvar_;vi++)
		for (size_t vi=0;vi<1;vi++) // Only calc x-component
		{
			for (size_t i=0;i<inmesh.nx_;i++)
			{
				for (size_t j=0;j<inmesh.ny_;j++)
				{
					for (size_t k=0;k<inmesh.nz_;k++)
					{
						energy += inmesh(i,j,k,vi)*inmesh(i,j,k,vi);
					}
				}
			}
		}
		//Real size = inmesh.nx_*inmesh.ny_*inmesh.nz_*inmesh.nvar_;
		Real size = inmesh.nx_*inmesh.ny_*inmesh.nz_;
		energy = energy/size;
	}
	
        void calc_mesh_umax(const Mesh &inmesh)
        {
                umax = 0.0;
		//for (size_t vi=0;vi<inmesh.nvar_;vi++)
		for (size_t vi=0;vi<1;vi++) // Only x-component
		{
			for (size_t i=0;i<inmesh.nx_;i++)
			{
				for (size_t j=0;j<inmesh.ny_;j++)
				{
					for (size_t k=0;k<inmesh.nz_;k++)
					{
						if (fabs(inmesh(i,j,k,vi)) >= umax)
							umax = fabs(inmesh(i,j,k,vi));
					}
				}
			}
		}
				
        }

        /// ---------------------------------------------------------
            
}; //End class Stats
