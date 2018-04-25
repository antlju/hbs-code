#include "diffops.h"

/// Computes omega = curl(u).
/// Where u is a 3-vector bundle and omega is a 3-vector pencil.
void curl(const Bundle &u, Pencil &omega, const Real xfac, const Real yfac, const Real zfac)
{
        Real dx,dy,dz;
        for (size_t i=0;i<u.nx_;i++)
        {
                
                /// Compute omega_1 = d_2u_3-d_3u_2
                dy = dely(u,yfac,i,2);
                dz = delz(u,zfac,i,1); 
                omega(i,0,0) = dy-dz;
                
                /// Compute omega_2 = d_3u_1-d_1u_3
                dz = delz(u,zfac,i,0);
                dx = delx(u,xfac,i,2);
                omega(i,0,1) = dz-dx;
                        
                /// Compute omega_3 = d_1u_2-d_2u_1
                dx = delx(u,xfac,i,1);
                dy = dely(u,yfac,i,0);
                omega(i,0,2) = dx-dy;
       
        }

}

