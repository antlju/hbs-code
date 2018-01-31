#include "common.h"
#include "utils.h"
#include "diffops.h"

void scalfieldsetup(Rfield U, const Real *x1, const Real *x2, const Real *x3)
{
        Real k=1.0,A=1.0,B=1.0,C=1.0;
        
        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                        {
                                U[vfidx(xi,yi,zi,0)] = sin(x2[yi-NGHOST]);
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
        Real dx = (L1-L0)/(NX-1);
        
        //Make 1d grids
        Real *x = linspace(NX,L0,L1);
        Real *z = linspace(NY,L0,L1);
        Real *y = linspace(NZ,L0,L1);
        
        Rfield u(1);

        scalfieldsetup(u,x,y,z);
        periodic_bc(u); 

        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                cout << u[vfidx(NGHOST,xi,NGHOST)] << "\t" << cos(x[xi-NGHOST]) << endl;
        }
        cout << endl;
        for (Int xi=0;xi<NX+NGHOST+2;xi++)
        {
                cout << u[vfidx(NGHOST,xi,NGHOST)] << endl;
        }
        cout << endl;
        Rfield g(3);
        grad(g,u,1.0/(dx),1.0/(dx),1.0/(dx));

        
        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                        {
                                cout << g[vfidx(NGHOST,yi,NGHOST,1)] << "\t";
                        }
                        cout << endl;
                }
                cout << endl;
        }
        
        
        return 0;
}
