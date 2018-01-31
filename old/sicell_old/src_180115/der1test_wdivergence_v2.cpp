#include "common.h"

Real *linspace(Int size, Real start, Real end)
{
        Real dx = (end-start)/size;
        Real *axis = new Real[size];
        for (Int i=0;i<size;i++)
        {
                axis[i] = start+dx*i;
                //cout << axis[i] << endl;
        }
        //cout << endl;
        return axis;
}

void setup(Rfield U, Real *xx,Real *yy,Real *zz)
{
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //set u_1 = sin z+cos y
                                U[vfidx(xi,yi,zi,0)] = 1.0;//sin(zz[zi-NGHOST])+cos(yy[yi-NGHOST]);

                                //set u_2 = sin x+cos z
                                U[vfidx(xi,yi,zi,1)] = sin(xx[xi-NGHOST])+cos(zz[zi-NGHOST]);

                                //set u_3 = sin y+cos x
                                U[vfidx(xi,yi,zi,2)] = sin(yy[yi-NGHOST])+cos(xx[xi-NGHOST]);
                                
                        }    
                }    
        }
}

void periodic_bc(Rfield U)
{
        for (Int vi=0;vi<U.nvar();vi++)
        {
                for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                {
                        for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                        {
                                for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                                {
                                        //x PBCs
                                        U[vfidx(0,yi,zi,vi)] = U[vfidx(NX+NGHOST-2,yi,zi,vi)];
                                        U[vfidx(1,yi,zi,vi)] = U[vfidx(NX+NGHOST-1,yi,zi,vi)];
                                        U[vfidx(NX+NGHOST,yi,zi,vi)] = U[vfidx(NGHOST+1,yi,zi,vi)];
                                        U[vfidx(NX+NGHOST+1,yi,zi,vi)] = U[vfidx(NGHOST+2,yi,zi,vi)];

                                        //y PBCs
                                        U[vfidx(xi,0,zi,vi)] = U[vfidx(xi,NY+NGHOST-2,zi,vi)];
                                        U[vfidx(xi,1,zi,vi)] = U[vfidx(xi,NY+NGHOST-1,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST,zi,vi)] = U[vfidx(xi,NGHOST+1,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST+1,zi,vi)] = U[vfidx(xi,NGHOST+2,zi,vi)];
                                
                                        //z PBCs
                                        U[vfidx(xi,yi,0,vi)] = U[vfidx(xi,yi,NZ+NGHOST-2,vi)];
                                        U[vfidx(xi,yi,1,vi)] = U[vfidx(xi,yi,NZ+NGHOST-1,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST,vi)] = U[vfidx(xi,yi,NGHOST+1,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST+1,vi)] = U[vfidx(xi,yi,NGHOST+1,vi)];
                              
                                }    
                        }    
                }
        }
}

int main()
{
        //Params
        //Real k=1.0,A=1.0,B=1.0,C=1.0;
        
        //Define grid
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/NX;
        
        //Make space for ABC field "u" and g.
        Rfield u(3);
        Rfield g(1);
        
        //Make 1d grid
        Real *x = linspace(NX,L0,L1);
        Real *z = linspace(NY,L0,L1);
        Real *y = linspace(NZ,L0,L1);

        setup(u,x,y,z);
        periodic_bc(u);

        //Try to take d_x (sin z+cos y) (should be zero for all)
        Real dxdu1;
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                dxdu1 = (1.0/dx) * der1(u[vfidx(xi-2,yi,zi,0)],u[vfidx(xi-1,yi,zi,0)],
                                                      u[vfidx(xi+1,yi,zi,0)],u[vfidx(xi+2,yi,zi,0)]);
                                //cout << u[vfidx(xi,yi,zi,0)] << "\t";
                                cout << dxdu1 << "\t";
                                
                        }
                        cout << endl;
                }
                cout << endl;
        }
}
