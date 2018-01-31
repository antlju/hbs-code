#include "common.h"
#include "diffops.h"

//u = grad(f), where u is a vector field and f is a scalar field.
void grad(Rfield u, const Rfield f, const Real xfactor, const Real yfactor, const Real zfactor)
{
        for (Int zi = 0;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=0;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=0;xi<NX+NGHOST;xi++)
                        {
                                //u_1 = d_x f
                                u[vfidx(xi,yi,zi,0)] = xfactor * der1(f[vfidx(xi-2,yi,zi)],f[vfidx(xi-1,yi,zi)],
                                                            f[vfidx(xi+1,yi,zi)],f[vfidx(xi+2,yi,zi)]);
                                
                                //u_2 = d_y f
                                u[vfidx(xi,yi,zi,1)] = yfactor * der1(f[vfidx(xi,yi-2,zi)],f[vfidx(xi,yi-1,zi)],
                                                            f[vfidx(xi,yi+1,zi)],f[vfidx(xi,yi+2,zi)]);

                                //u_3 = d_z f
                                u[vfidx(xi,yi,zi,2)] = zfactor * der1(f[vfidx(xi,yi,zi-2)],f[vfidx(xi,yi,zi-1)],
                                                            f[vfidx(xi,yi,zi+1)],f[vfidx(xi,yi,zi+2)]);
                        }
                }
        }
}
