#include "common.h"
#include "utils.h"
#include "diffops.h"

void fieldsetup(Rfield uABC, const Real *x1, const Real *x2, const Real *x3)
{
        Real k=3.0,A=1.0,B=1.0,C=1.0;
        
        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                        {
                                //set u_1(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,0)] = exp(k*x1[xi-NGHOST]); //exp(k*x1[xi-NGHOST]);
                                        //A*sin(k*x3[zi-NGHOST])+C*cos(k*x2[yi-NGHOST]);

                                //set u_2(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,1)] = exp(k*x2[yi-NGHOST]);
                                        //B*sin(k*x1[xi-NGHOST])+A*cos(k*x3[zi-NGHOST]);

                                //set u_3(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,2)] = exp(k*x3[zi-NGHOST]);
                                        //C*sin(k*x2[yi-NGHOST])+B*cos(k*x1[xi-NGHOST]);
                        }
                }
        }
}

Int main()
{
//Params
        //Real k=1.0;//,A=1.0,B=1.0,C=1.0;
        
        //Define grid
        Real L0 = 0.0;
        Real L1 = 1.0;
        Real dx = (L1-L0)/(NX);
        
        //Make 1d grids
        Real *x = linspace(NX,L0,L1);
        Real *z = linspace(NY,L0,L1);
        Real *y = linspace(NZ,L0,L1);
        
        Rfield u(3);

        fieldsetup(u,x,y,z);
        periodic_bc(u); 

        Rfield g(1);
        div(g,u,1.0/(dx),1.0/(dx),1.0/(dx));

        Int yi = NGHOST;
        Int zi = NGHOST;
        Int xi = NGHOST;
        Real k = 1.0;
/*
  for (xi=NGHOST;xi<NX+NGHOST;xi++)
  {
  Real d1u1 = (1.0/dx) *
  der1(u[vfidx(xi-2,yi,zi,0)],u[vfidx(xi-1,yi,zi,0)],
  u[vfidx(xi+1,yi,zi,0)],u[vfidx(xi+2,yi,zi,0)]);
  cout << u[vfidx(xi,yi,zi,0)] << "\t" << d1u1 << "\t" << g[vfidx(xi,yi,zi)] << endl;
  //cout << u[vfidx(xi,yi,zi,0)] << endl;
                
  }
*/
        /*
        for (yi=NGHOST;yi<NY+NGHOST;yi++)
        {
                Real d2u2 = (1.0/dx) *
                        der1(u[vfidx(xi,yi-2,zi,1)],u[vfidx(xi,yi-1,zi,1)],
                             u[vfidx(xi,yi+1,zi,1)],u[vfidx(xi,yi+2,zi,1)]);
                cout << u[vfidx(xi,yi,zi,1)] << "\t" << d2u2 << "\t" << g[vfidx(xi,yi,zi)] << endl;
                //cout << u[vfidx(xi,yi,zi,0)] << endl;
                
        }
        */

        
        for (xi=NGHOST+2;xi<NX+NGHOST-2;xi++)
        {
                for (yi=NGHOST+2;yi<NY+NGHOST-2;yi++)
                {
                        for (zi=NGHOST+2;zi<NZ+NGHOST-2;zi++)
                        {
                                /*
                                cout << g[vfidx(xi,yi,zi)] << "\t"
                                     << k*(u[vfidx(xi,yi,zi,0)]
                                           +u[vfidx(xi,yi,zi,1)]
                                           +u[vfidx(xi,yi,zi,2)]) << endl;
                                */
                                cout << g[vfidx(xi,yi,zi)]/(u[vfidx(xi,yi,zi,0)]
                                           +u[vfidx(xi,yi,zi,1)]
                                                            +u[vfidx(xi,yi,zi,2)]) << endl;
                                
                        }
                }
        }
                 

        /*
        cout << u[vfidx(NGHOST+2,NGHOST+2,NGHOST+2,0)] << endl;
        cout << u[vfidx(NGHOST+2,NGHOST+2,NGHOST+2,1)] << endl;
        cout << u[vfidx(NGHOST+2,NGHOST+2,NGHOST+2,2)] << endl;
        cout << g[vfidx(NGHOST+2,NGHOST+2,NGHOST+2)] << endl;
        */
        
        return 0;
}
