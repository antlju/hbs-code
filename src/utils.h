#pragma once
#include "common.h"

void printfield(const Field4 &u) //Prints contents to command line (excluding ghost points)
{
        for (size_t vii=0;vii<u.nvar_;vii++)
        {
                std::cout << "------- component: " << vii << " -------" << std::endl;
                std::cout << "---------------------------------------" << std::endl;
                for (size_t ii=0;ii<u.nx_;ii++)
                {
                        for (size_t jj=0;jj<u.ny_;jj++)
                        {
                                for (size_t kk=0;kk<u.nz_;kk++)
                                {
                                        //std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii,jj,kk,vii) << std::endl;
                                        std::cout << u(ii,jj,kk,vii) << "\t";
                                }
                                std::cout << "\n";
                        }
                        std::cout << "---------------------------------------" << std::endl;
                        std::cout << std::endl;
                }
                        
        }
}

void linspace(sPencil &u, Real start, Real end, Real dx)
{
        Real size = u.nx_;
        //Real dx = (end-start)/size;
        for (size_t i=0;i<size;i++)
        {
                u(i) = start+dx*i;
                //std::cout << u(i) << std::endl;
        }
        //std::cout << std::endl;
}

Field4 pwiseffdiv(const Field4 &u1,const Field4 &u2)
{
        Field4 div(u1.nx_,u1.ny_,u1.nz_,u1.nvar_);
        for (size_t vii=0;vii<u1.nvar_;vii++)
        {
                for (size_t ii=0;ii<u1.nx_;ii++)
                {
                        for (size_t jj=0;jj<u1.ny_;jj++)
                        {
                                for (size_t kk=0;kk<u1.nz_;kk++)
                                {
                                        div(ii,jj,kk,vii) = u1(ii,jj,kk,vii)/u2(ii,jj,kk,vii);   
                                }
              
                        }
              
                }
                        
        }

        return div;     
}

void apply_pbc(Field4 &u);
