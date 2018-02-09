#include "common.h"

void apply_pbc(Field4 &u)
{
        size_t Ng = u.ng_,Nvars=u.nvar_;
        size_t Nx=u.nx_,Ny=u.ny_,Nz=u.nz_;
        
        for (size_t vi=0;vi<Nvars;vi++)
        {
                
                for (size_t i=0;i<Nx;i++)
                {
                        for (size_t j=0;j<Ny;j++)
                        {
                                for (size_t k=0;k<Nz;k++)
                                {
                                        for (size_t l=0;l<Ng;l++)
                                        {
                                                //set pbc along x
                                                u(l-Ng,j,k,vi) = u(Nx-(Ng-l),j,k,vi);
                                                u(Nx+l,j,k,vi) = u(l,j,k,vi);

                                                //set pbc along y
                                                u(i,l-Ng,k,vi) = u(i,Ny-(Ng-l),k,vi);
                                                u(i,Ny+l,k,vi) = u(i,l,k,vi);

                                                //set pbc along z
                                                u(i,j,l-Ng,vi) = u(i,j,Nz-(Ng-l),vi);
                                                u(i,j,Nz+l,vi) = u(i,j,l,vi);
                                        }
                                }
                        }
                }
               
        }
}
