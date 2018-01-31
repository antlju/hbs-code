#include "common.h"
#include "diffops.h"

//omega = curl(U). 
void curl(Rfield omega, const Rfield U, const Real xfactor, const Real yfactor, const Real zfactor)
{
        Real d1,d2;
        
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //Compute omega_1 = d_2u_3-d_3u_2
                                d1 = yfactor * der1(U[vfidx(xi,yi-2,zi,2)],U[vfidx(xi,yi-1,zi,2)],
                                                    U[vfidx(xi,yi+1,zi,2)],U[vfidx(xi,yi+2,zi,2)]);
                                
                                d2 = zfactor * der1(U[vfidx(xi,yi,zi-2,1)],U[vfidx(xi,yi,zi-1,1)],
                                                    U[vfidx(xi,yi,zi+1,1)],U[vfidx(xi,yi,zi+2,1)]);
                                
                                omega[vfidx(xi,yi,zi,0)] = d1-d2;
                                      
                                //Compute omega_2 = d_3u_1-d_1u_3
                                d1 = zfactor * der1(U[vfidx(xi,yi,zi-2,0)],U[vfidx(xi,yi,zi-1,0)],
                                                    U[vfidx(xi,yi,zi+1,0)],U[vfidx(xi,yi,zi+2,0)]);
                                
                                d2 = xfactor * der1(U[vfidx(xi-2,yi,zi,2)],U[vfidx(xi-1,yi,zi,2)],
                                                    U[vfidx(xi+1,yi,zi,2)],U[vfidx(xi+2,yi,zi,2)]);

                                omega[vfidx(xi,yi,zi,1)] = d1-d2;
                                
                                //Compute omega_3 = d_1u_2-d_2u_1
                                d1 = xfactor * der1(U[vfidx(xi-2,yi,zi,1)],U[vfidx(xi-1,yi,zi,1)],
                                                    U[vfidx(xi+1,yi,zi,1)],U[vfidx(xi+2,yi,zi,1)]);

                                d2 = yfactor * der1(U[vfidx(xi,yi-2,zi,0)],U[vfidx(xi,yi-1,zi,0)],
                                                    U[vfidx(xi,yi+1,zi,0)],U[vfidx(xi,yi+2,zi,0)]);

                                omega[vfidx(xi,yi,zi,2)] = d1-d2;
                                
                        }
                }
        }

        
}
