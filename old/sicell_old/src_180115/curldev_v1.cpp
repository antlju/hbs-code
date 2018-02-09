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
        for (int xi=NGHOST;xi<NX+NGHOST;xi++)
        {
                for (int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
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

void grad(Rfield U,Int xi,Int yi,Int zi,Real *uij,Real factor)
{
        

        for (Int i=0;i<3;i++)
        {
                uij[mat2lin(i,0)] = factor * der1(U[vfidx(xi-2,yi,zi,i)],U[vfidx(xi-1,yi,zi,i)],
                                                  U[vfidx(xi+1,yi,zi,i)],U[vfidx(xi+2,yi,zi,i)]);
                        
                uij[mat2lin(i,1)] = factor * der1(U[vfidx(xi,yi-2,zi,i)],U[vfidx(xi,yi-1,zi,i)],
                                                  U[vfidx(xi,yi+1,zi,i)],U[vfidx(xi,yi+1,zi,i)]);
                        
                uij[mat2lin(i,2)] = factor * der1(U[vfidx(xi,yi,zi-2,i)],U[vfidx(xi,yi,zi-1,i)],
                                                  U[vfidx(xi,yi,zi+1,i)],U[vfidx(xi,yi,zi+2,i)]);
        }
        
}


//omega = curl(U). 
void curl(Rfield omega, const Rfield U, const Real xfactor)
{

        Real Uij[9];
        
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                                grad(U,xi,yi,zi,Uij,xfactor);
                                omega[vfidx(xi,yi,zi,0)] = Uij[mat2lin(2,1)]-Uij[mat2lin(1,2)]; //omega_1 = d_2u_3-d_3u_2
                                omega[vfidx(xi,yi,zi,1)] = Uij[mat2lin(0,2)]-Uij[mat2lin(2,0)];
                                omega[vfidx(xi,yi,zi,2)] = Uij[mat2lin(1,0)]-Uij[mat2lin(0,1)];
                        }
                }
        }

        
}
/*
void div(Rfield g, const Rfield U, const Real xfactor)
{

        Real Uij[3][3];
        
        for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
        {
                for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                {
                        for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                        {
                               
                                grad(U,xi,yi,zi,Uij,xfactor);
                                g[vfidx(xi,yi,zi)] = Uij[0][0]+Uij[1][1]+Uij[2][2];

                        }
                }
        }

        
}
*/
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
                                        U[vfidx(NX+NGHOST,yi,zi,vi)] = U[vfidx(NGHOST,yi,zi,vi)];
                                        U[vfidx(NX+NGHOST+1,yi,zi,vi)] = U[vfidx(NGHOST+1,yi,zi,vi)];

                                        //y PBCs
                                        U[vfidx(xi,0,zi,vi)] = U[vfidx(xi,NY+NGHOST-2,zi,vi)];
                                        U[vfidx(xi,1,zi,vi)] = U[vfidx(xi,NY+NGHOST-1,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST,zi,vi)] = U[vfidx(xi,NGHOST,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST+1,zi,vi)] = U[vfidx(xi,NGHOST+1,zi,vi)];
                                
                                        //z PBCs
                                        U[vfidx(xi,yi,0,vi)] = U[vfidx(xi,yi,NZ+NGHOST-2,vi)];
                                        U[vfidx(xi,yi,1,vi)] = U[vfidx(xi,yi,NZ+NGHOST-1,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST,vi)] = U[vfidx(xi,yi,NGHOST,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST+1,vi)] = U[vfidx(xi,yi,NGHOST+1,vi)];
                              
                                }    
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
        
        //Make space for ABC field "u" and define components.
        Rfield u(3);
        Rfield curluAn(3);
        //Rfield ux = u.subfield(0,1);
        //Rfield uy = u.subfield(1,1);
        //Rfield uz = u.subfield(2,1);

        //Make 1d grid
        Real *x = linspace(NX,L0,L1);
        Real *z = linspace(NY,L0,L1);
        Real *y = linspace(NZ,L0,L1);

        ABC_setup(u,k,A,B,C,x,y,z);
        periodic_bc(u);
        //curlABCanalytic(curluAn,k,A,B,C,x,y,z);

        //Make space for curl(ABC)
        Rfield cu(3);
        //Rfield cux = cu.subfield(0,1);
        //Rfield cuy = cu.subfield(1,1);
        //Rfield cuz = cu.subfield(2,1);
        Rfield g(1);
        
        curl(cu,u,1.0/dx);
        //div(g,u,1.0/dx);
        
        Rfield udotcu(1);
        Rfield udotu(1);
        dotprod(u,cu,udotcu);
        dotprod(u,u,udotu);

        
        //Print some Rfield
        for (Int vi=0;vi<1;vi++)
        {
                for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                {
                        cout << "\tzi = " << zi-NGHOST << "\t z = " << z[zi-NGHOST] << endl;
                                                        cout << "---------" << endl;
                        for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                        {
                                cout << "\t y = " << y[yi-NGHOST] << endl;
                                for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                                {
                                        //cout << "\t" << fabs(cu[vfidx(xi,yi,zi,vi)]-curluAn[vfidx(xi,yi,zi,vi)]) << "\t";
                                        //cout << "\t" << curluAn[vfidx(xi,yi,zi,vi)] << "\t" << cu[vfidx(xi,yi,zi,vi)];
                                        cout.precision(17);
                                        //cout << "\t x = " << x[xi-NGHOST];
                                        if (abs(udotu[vfidx(xi,yi,zi,0)]) < 0.00001)
                                        {
                                                cout << "\t" << "1";
                                        }
                                        else
                                        {
                                                cout << "\t" << udotu[vfidx(xi,yi,zi,0)];
                                                
                                        }
                                        
                                        //cout << "\t" << g[vfidx(xi,yi,zi)];
                                }
                                cout << endl;
                        }
                        cout << endl;
                }
                cout << "------------ end component --------------- " << endl;
        }
        
        
        return 0;
}
