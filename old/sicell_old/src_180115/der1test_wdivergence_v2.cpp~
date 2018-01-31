#include "common.h"

Real *linspace(Int size, Real start, Real end)
{
        Real dx = (end-start)/size;
        Real *axis = new Real[size];
        for (Int i=0;i<size;i++)
        {
                axis[i] = start+dx*i;
                cout << axis[i] << endl;
        }
        cout << endl;
        return axis;
}


void ABC_setup(Rfield uABC, const Real k, const Real A, const Real B, const Real C, const Real *x1, const Real *x2, const Real *x3)
{
        for (int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //set u_1(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,0)] = 
                                        A*sin(k*x3[zi-NGHOST])+C*cos(k*x2[yi-NGHOST]);

                                //set u_2(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,1)] =
                                        B*sin(k*x1[xi-NGHOST])+A*cos(k*x3[zi-NGHOST]);

                                //set u_3(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,2)] =
                                        C*sin(k*x2[yi-NGHOST])+B*cos(k*x1[xi-NGHOST]);
                        }
                }
        }
        
}

void curlABCanalytic(Rfield uABC, const Real k, const Real A, const Real B, const Real C, const Real *x1, const Real *x2, const Real *x3)
{
        for (int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //set u_1(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,0)] = 
                                        A*k*cos(k*x2[yi-NGHOST])-(-C*k*sin(k*x3[zi-NGHOST]));

                                //set u_2(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,1)] =
                                        B*k*cos(k*x3[yi-NGHOST])-(-A*k*sin(k*x1[xi-NGHOST]));

                                //set u_3(xi,yi,zi)
                                uABC[vfidx(xi,yi,zi,2)] =
                                        C*k*cos(k*x1[xi-NGHOST])-(-B*k*sin(k*x2[yi-NGHOST]));
                        }
                }
        }
        
}

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
                                //Compute omega_1 = d_2u_3-d_3 u_2
                                d1 = yfactor * der1(U[vfidx(xi,yi-2,zi,2)],U[vfidx(xi,yi-1,zi,2)],
                                                    U[vfidx(xi,yi+1,zi,2)],U[vfidx(xi,yi+2,zi,2)]);
                                
                                d2 = zfactor * der1(U[vfidx(xi,yi,zi-2,1)],U[vfidx(xi,yi,zi-1,1)],
                                                    U[vfidx(xi,yi,zi+1,1)],U[vfidx(xi,yi,zi+2,1)]);
                                
                                omega[vfidx(xi,yi,zi,0)] = d1-d2;
                                      
                                //Compute omega_2 = d_3u_1-d_1u_3
                                d1 = xfactor * der1(U[vfidx(xi,yi,zi-2,0)],U[vfidx(xi,yi,zi-1,0)],
                                                    U[vfidx(xi,yi,zi+1,0)],U[vfidx(xi,yi,zi+2,0)]);
                                
                                d2 = zfactor * der1(U[vfidx(xi-2,yi,zi,2)],U[vfidx(xi-1,yi,zi,2)],
                                                    U[vfidx(xi+1,yi,zi,2)],U[vfidx(xi+2,yi,zi,2)]);

                                omega[vfidx(xi,yi,zi,1)] = d1-d2;
                                
                                //Compute omega_3 = d_1u_2-d_2u_1
                                d1 = yfactor * der1(U[vfidx(xi-2,yi,zi,1)],U[vfidx(xi-1,yi,zi,1)],
                                                    U[vfidx(xi+1,yi,zi,1)],U[vfidx(xi+2,yi,zi,1)]);

                                d2 = xfactor * der1(U[vfidx(xi,yi-2,zi,0)],U[vfidx(xi,yi-1,zi,0)],
                                                    U[vfidx(xi,yi+1,zi,0)],U[vfidx(xi,yi+2,zi,0)]);

                                omega[vfidx(xi,yi,zi,2)] = d1-d2;
                                
                        }
                }
        }

        
}

void dotprod(Rfield a, Rfield b, Rfield c)
{

        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                c[vfidx(xi,yi,zi,0)] =
                                        a[vfidx(xi,yi,zi,0)]*b[vfidx(xi,yi,zi,0)]
                                        +a[vfidx(xi,yi,zi,1)]*b[vfidx(xi,yi,zi,1)]
                                        +a[vfidx(xi,yi,zi,2)]*b[vfidx(xi,yi,zi,2)];
                                        
                        }
                }
        }

}

//g = div(a),
void div(Rfield u, Rfield g, Real xfactor, Real yfactor, Real zfactor)
{
        Real dxdu1,dydu2,dzdu3;
        
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                //compute d_x u_1 at the point (xi,yi,zi)
                                dxdu1 = xfactor * der1(u[vfidx(xi-2,yi,zi,0)],u[vfidx(xi-1,yi,zi,0)],
                                             u[vfidx(xi+1,yi,zi,0)],u[vfidx(xi+2,yi,zi,0)]);
                                
                                //compute d_y u_2 at the point (xi,yi,zi)
                                dydu2 = yfactor * der1(u[vfidx(xi,yi-2,zi,1)],u[vfidx(xi,yi-1,zi,1)],
                                             u[vfidx(xi,yi+1,zi,1)],u[vfidx(xi,yi+2,zi,1)]);

                                //compute d_z u_3 at the point (xi,yi,zi)
                                dzdu3 = zfactor * der1(u[vfidx(xi,yi,zi-2,2)],u[vfidx(xi,yi,zi-1,2)],
                                             u[vfidx(xi,yi,zi+1,2)],u[vfidx(xi,yi,zi+2,2)]);

                                //set g(xi,yi,zi)
                                g[vfidx(xi,yi,zi,0)] = dxdu1+dydu2+dzdu3;
                        }
                }
        }
}

int main()
{
        //Params
        Real k=1.0,A=1.0,B=1.0,C=1.0;
        
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

        ABC_setup(u,k,A,B,C,x,y,z);

        div(u,g,1.0/dx,1.0/dx,1.0/dx);
        
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                cout << g[vfidx(xi,yi,zi,0)] << "\t";
                        }
                        cout << endl;
                }
                cout << endl;
        }
        return 0;
}
