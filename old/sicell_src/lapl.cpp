#include "common.h"
#include "diffops.h"

//d2f = Laplacian(f).

void lapl(Rfield d2f, const Rfield f, const Real xfactor, const Real yfactor, const Real zfactor)
{
        for (Int vi = 0;vi<d2f.nvar();vi++)
        {
                for (Int zi = 0;zi<NZ+NGHOST;zi++)
                {
                        for (Int yi=0;yi<NY+NGHOST;yi++)
                        {
                                for (Int xi=0;xi<NX+NGHOST;xi++)
                                {
                                        //We sum up the del^2_i f in this variable and store it at the appropriate point in the scalar field d2f.
                                        Real valAtPoint = xfactor
                                                * der2(f[vfidx(xi-2,yi,zi,vi)],f[vfidx(xi-1,yi,zi,vi)],
                                                       f[vfidx(xi,yi,zi,vi)],
                                                       f[vfidx(xi+1,yi,zi,vi)],f[vfidx(xi+2,yi,zi,vi)]);

                                        valAtPoint += yfactor
                                                * der2(f[vfidx(xi,yi-2,zi,vi)],f[vfidx(xi,yi-1,zi,vi)],
                                                       f[vfidx(xi,yi,zi,vi)],
                                                       f[vfidx(xi,yi+1,zi,vi)],f[vfidx(xi,yi+2,zi,vi)]);

                                        valAtPoint += zfactor
                                                * der2(f[vfidx(xi,yi,zi-2,vi)],f[vfidx(xi,yi,zi-1,vi)],
                                                       f[vfidx(xi,yi,zi,vi)],
                                                       f[vfidx(xi,yi,zi+1,vi)],f[vfidx(xi,yi,zi+2,vi)]);

                                        d2f[vfidx(xi,yi,zi,vi)] = valAtPoint;
                                }
                        }
                }
        }
}
