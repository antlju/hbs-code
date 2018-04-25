#include "includes.h"

void printfield(const Mesh &u) //Prints contents to command line (excluding ghost points)
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
