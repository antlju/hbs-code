#include "includes.h"

#define DD 1.0 // l^2/s
#define MAXTSTEPS 200

std::string set_fname(std::string base, std::string end, Int step)
{
        std::string timeval = "_t_";
        std::string timeNo = "_tNum_";
        return base+timeNo+std::to_string(step)+end;
}

Int writeToFile(std::string fname, const Mesh &u, const Pencil &x, const Pencil &y, const Pencil &z)
{
        std::ofstream openfile("/home/anton/dev/hbs/simdata/"+fname, std::ios::trunc);
        openfile << "i,j,k,x,y,z,f" << std::endl;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                openfile << i << "," << j << "," << k << ","
                                         << x(i) << "," << y(j) << "," << z(k) << ","
                                         << u(i,j,k,0) << std::endl;
                        }
                        
                }
        }
        openfile.close();
        std::cout << "Wrote " << fname << std::endl;
        return 0;

}

void gauss_init(Mesh &u,const Pencil &x, const Pencil &y, const Pencil &z,const Real dx)
{
        Real xshift=(((u.nx_/2))*dx),yshift=((u.ny_/2)*dx),zshift=((u.nz_/2)*dx);
        Real std2 = pow(0.7,2);
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                Real xx = x(i)-xshift;
                                Real yy = y(j)-yshift;
                                Real zz = z(i)-zshift;
                                Real r2 = pow(xx,2)+pow(yy,2)+pow(zz,2);
                                
                                Real gaussian = exp(-0.5*(r2/std2));
                                //std::cout << gaussian << std::endl;
                                u(i,j,k,0) = gaussian;
                                
                                /*
                                u(u.nx_/2+i,u.ny_/2+j,u.nz_/2+k) = 0.8;
                                
                                u(u.nx_/2-i,u.ny_/2+j,u.nz_/2+k) = 0.8;
                                u(u.nx_/2+i,u.ny_/2-j,u.nz_/2+k) = 0.8;
                                u(u.nx_/2+i,u.ny_/2+j,u.nz_/2-k) = 0.8;
                                
                                u(u.nx_/2-i,u.ny_/2-j,u.nz_/2-k) = 0.8;
                                
                                u(u.nx_/2+i,u.ny_/2-j,u.nz_/2-k) = 0.8;
                                u(u.nx_/2-i,u.ny_/2+j,u.nz_/2-k) = 0.8;
                                u(u.nx_/2-i,u.ny_/2-j,u.nz_/2+k) = 0.8;
                                */
                        }
                } 
        }
        //u(u.nx_/2,u.ny_/2,u.nz_/2,0) = 1.0;
}


void diffops(Mesh &ff, Mesh &dff, const Real xfac, const Real yfac, const Real zfac)
{
        Bundle LaplB(ff.nx_,1); /// scalar bundle
        Pencil LaplP(ff.nx_,1); /// scalar pencil
        
        for (size_t j=0;j<ff.ny_;j++)
        {
                for (size_t k=0;k<ff.nz_;k++)
                {
                        ff2bundle(ff,LaplB,j,k);
                        lapl(LaplB,LaplP,xfac,yfac,zfac);
                        pencil2ff(LaplP,dff,j,k);
                }
        }
}

void RK2(Mesh &ff, Mesh &dff, const Real dt, const Int maxtsteps, const Real dx, const Real dy, const Real dz,const Pencil &x)
{
        Mesh k1(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k2(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh ff1(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);

        for (Int tstep=1;tstep<maxtsteps;tstep++)
        {
                diffops(ff,dff,dx,dy,dz); /// Calculate k1 = dt*f(t_n,y_n)
                //dff.print();
                k1 = dff*dt;
                //k1.print();
                
                ff1 = ff+(k1*0.5); /// Calculate y_n + 1/2*k1
                //ff1.print();
                diffops(ff1,dff,dx,dy,dz); //Calculate k2 = dt*f(t_n+0.5*dt,y_n+0.5*k1)
                k2 = dff*dt;
                //k2.print();
        
                ff = ff + k2; /// Step with error of O(dt^3)
                
                //ff.print();
                writeToFile(set_fname("rk2test_gaussInit_2",".csv",tstep),ff,x,x,x);
        }
        //std::cout << "completed timestep #" << tstep << std::endl;
        
     
        
}

void RK4(Mesh &ff, Mesh &dff, const Real dt, const Int maxtsteps, const Real dx, const Real dy, const Real dz,const Pencil &x)
{
        Mesh k1(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k2(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k3(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh k4(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);
        Mesh ff1(ff.nx_,ff.ny_,ff.nz_,ff.nvar_);

        Real c1=1.0/6.0,c2=1.0/3.0;

                for (Int tstep=1;tstep<maxtsteps;tstep++)
        {
                diffops(ff,dff,dx,dy,dz);
                k1 = dff*dt;

                ff1 = ff+(k1*0.5);
                diffops(ff1,dff,dx,dy,dz);
                k2 = dff*dt;

                ff1 = ff+(k2*0.5);
                diffops(ff1,dff,dx,dy,dz);
                k3 = dff*dt;

                ff1 = ff+(k3*0.5);
                diffops(ff1,dff,dx,dy,dz);
                k4 = dff*dt;

                ff = ff + k1*c1 + k2*c2 + k3*c2 + k4*c1;
                writeToFile(set_fname("rk4test_gaussInit",".csv",tstep),ff,x,x,x);
        }
        //std::cout << "completed timestep #" << tstep << std::endl;
        
     
        
}


Int main()
{
        Int Nx=NX,Ny=NY,Nz=NZ;
        Int maxtstep = MAXTSTEPS;
        
        Real L0=0;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/Nx;
        Real invdx2 = 1.0/pow(dx,2);

        Pencil x(Nx);
        linspace(x,L0,L1,dx);
        
        Real diffC = DD;
        Real dt = 0.1*pow(dx,2)/diffC;

        Mesh ff(Nx,Ny,Nz,1);
        Mesh dff(Nx,Ny,Nz,1);
        gauss_init(ff,x,x,x,dx);
        //ff.print();

        //writeToFile(set_fname("rk2test_gaussInit_2",".csv",0),ff,x,x,x);
        //RK2(ff,dff,dt,maxtstep,invdx2,invdx2,invdx2,x);

        writeToFile(set_fname("rk4test_gaussInit",".csv",0),ff,x,x,x);
        RK4(ff,dff,dt,maxtstep,invdx2,invdx2,invdx2,x);
        return 0;
}
