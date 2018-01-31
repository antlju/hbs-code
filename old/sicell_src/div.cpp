#include "common.h"
#include "diffops.h"

//g = div(u)
//g is a scalar field and u is a vector field.

void div(Rfield g, const Rfield u, const Real xfactor, const Real yfactor, const Real zfactor)
{
        Real valatpt=0.0;
        Real d1u1=0.0,d2u2=0.0,d3u3=0.0;
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //We sum up the x,y,z-derivatives in valAtPoint
                                //and store it at the appropriate point in the scalar field g.
                                //compute d_z u_3 at the point (xi,yi,zi)
                                d3u3 = zfactor *
                                        der1(u[vfidx(xi,yi,zi-2,2)],u[vfidx(xi,yi,zi-1,2)],
                                             u[vfidx(xi,yi,zi+1,2)],u[vfidx(xi,yi,zi+2,2)]);

                                
                                //compute d_x u_1 at the point (xi,yi,zi)
                                d1u1 = xfactor *
                                        der1(u[vfidx(xi-2,yi,zi,0)],u[vfidx(xi-1,yi,zi,0)],
                                             u[vfidx(xi+1,yi,zi,0)],u[vfidx(xi+2,yi,zi,0)]);
                                
                                //compute d_y u_2 at the point (xi,yi,zi)
                                d2u2 = yfactor *
                                        der1(u[vfidx(xi,yi-2,zi,1)],u[vfidx(xi,yi-1,zi,1)],
                                             u[vfidx(xi,yi+1,zi,1)],u[vfidx(xi,yi+2,zi,1)]);
                                
                                //set g(xi,yi,zi)
                                
                                //cout << d1u1 << "\t" << u[vfidx(xi,yi,zi,0)] << endl;
                                //cout << d2u2 << "\t" << u[vfidx(xi,yi,zi,1)] << endl;
                                //cout << d3u3 << "\t" << u[vfidx(xi,yi,zi,2)] << endl;
                                //cout << endl;

                                /*
                                cout << "[" << d1u1+d2u2+d3u3 << "\t" << u[vfidx(xi,yi,zi,0)]+u[vfidx(xi,yi,zi,1)]+u[vfidx(xi,yi,zi,2)] << "]" << "\t";
                                */
                                valatpt = d1u1+d2u2+d3u3;
                                //cout << valatpt << endl;
                                g[vfidx(xi,yi,zi)] = valatpt;
                        }
                        //cout << endl;
                }
                //cout << endl;
        }
        
}
