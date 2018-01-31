#include "common.h"
#include "utils.h"
#include "diffops.h"

void fieldsetup(Rfield uABC, const Real *x1, const Real *x2, const Real *x3)
{
        Real k=2.0,A=1.0,B=1.0,C=1.0;
        
        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                        {
                                //set u_1(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,0)] = sin(k*x1[xi-NGHOST]);
                                        //A*sin(k*x3[zi-NGHOST])+C*cos(k*x2[yi-NGHOST]);

                                //set u_2(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,1)] = sin(k*x2[yi-NGHOST]);
                                        //B*sin(k*x1[xi-NGHOST])+A*cos(k*x3[zi-NGHOST]);

                                //set u_3(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,2)] = sin(k*x3[zi-NGHOST]);
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
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/(NX);
        
        //Make 1d grids
        Real *x = linspace(NX,L0,L1);
        Real *z = linspace(NY,L0,L1);
        Real *y = linspace(NZ,L0,L1);
        
        Rfield u(3);

        fieldsetup(u,x,y,z);
        periodic_bc(u); 

        Rfield g(3);
        lapl(g,u,1.0/(dx*dx),1.0/(dx*dx),1.0/(dx*dx));


        for (Int vi=0;vi<u.nvar();vi++)
        {
                for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                {
                        for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                        {
                                for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                                {
                                        if (u[vfidx(xi,yi,zi,vi)] < 1e-15){
                                                cout << "uzero" << "\t";
                                        }
                                        else
                                        {
                                                cout << g[vfidx(xi,yi,zi,vi)]/u[vfidx(xi,yi,zi,vi)] << "\t";
                                        }
                                }
                                cout << endl;
                        }
                        cout << endl;
                }
        }
        
        
        return 0;
}
