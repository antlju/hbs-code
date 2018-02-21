

#include "includes.h"

#define DD 1.0 // l^2/s
#define MAXTSTEPS 2000

void u_equals_u_plus_dtxdu(Mesh &u, Mesh &du, const Real dt)
{
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                u(i,j,k,0) = u(i,j,k,0)+dt*du(i,j,k,0);
                        }
                }
        }
};

void timestep(Mesh &u, Mesh &du, const Real dt, const Real xfac, const Real yfac, const Real zfac)
{
        Bundle LaplB(u.nx_,1); /// scalar bundle
        Pencil LaplP(u.nx_,1); /// scalar pencil
        
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {
                        ff2bundle(u,LaplB,j,k);
                        lapl(LaplB,LaplP,xfac,yfac,zfac);
                        pencil2ff(LaplP,du,j,k);
                }
        }
        u_equals_u_plus_dtxdu(u,du,dt);
        
}

Real initGauss(

void pointsrc_init(Mesh &u)
{
        u(u.nx_/2,u.ny_/2,u.nz_/2,0) = 1.0;
}

std::string set_fname(std::string base, std::string end, Int step)
{
        std::string timeval = "_t_";
        std::string timeNo = "_tNum_";
        return base+timeNo+std::to_string(step)+end;
}

Int writeToFile(std::string fname, const Mesh &u)
{
        std::ofstream openfile("/home/anton/dev/hbs/simdata/"+fname, std::ios::trunc);
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                openfile << i << "," << j << "," << k << "," << u(i,j,k,0) << std::endl;
                        }
                        
                }
        }
        openfile.close();
        std::cout << "Wrote " << fname << std::endl;
        return 0;

}


Int main()
{
        Int Nx=NX,Ny=NY,Nz=NZ;
        Int maxtstep = MAXTSTEPS;
        
        Real L0=0;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/Nx;
        Real invdx2 = 1.0/pow(dx,2);
        
        Real diffC = DD;
        Real dt = 0.01*pow(dx,2)/diffC;

        Mesh u(Nx,Ny,Nz,1);
        pointsrc_init(u);

        Mesh du(Nx,Ny,Nz,1);
        
        writeToFile(set_fname("u",".dat",0),u);
                    
        for (Int tstep=1;tstep<maxtstep;tstep++)
        {
                timestep(u,du,dt,invdx2,invdx2,invdx2);
                writeToFile(set_fname("u",".csv",tstep),u);
        }
        

        return 0;
}
              
