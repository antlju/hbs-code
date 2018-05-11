#include "includes.h"

/// Sets filename for use in writeToFile()
std::string set_fname(const std::string base, const std::string end, const Int step)
{
        std::string timeval = "_t_";
        std::string timeNo = "_tNum_";
        return base+timeNo+std::to_string(step)+end;
}

/// Function to write a 3D scalar field and corresponding (x,y,z)-coordinates to file.
Int writeToFile(std::string fname, const Mesh &u, const Int vi, const Pencil &x, const Pencil &y, const Pencil &z)
{
        std::ofstream openfile("/home/anton/dev/hbs/simdata/"+fname, std::ios::trunc);
        //openfile << "i,j,k,x,y,z,f" << std::endl;
        //for (size_t i=0;i<u.nx_;i++)
        //{
        openfile << "#---------------------------------" << std::endl;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                /*
                                  openfile << i << "," << j << "," << k << ","
                                  << x(i) << "," << y(j) << "," << z(k) << ","
                                  << u(i,j,k,vi) << std::endl;
                                */

                                openfile << u(i,j,k,vi) << "\t";
                        }
                        openfile << std::endl;       
                }
                openfile << "#---------------------------------" << std::endl;
        }
                // }
        openfile.close();
        std::cout << "Wrote " << fname << std::endl;
        return 0;

}

Int writeToFile_1DArr(const std::string fname, const Mesh &u, const Int vi, const Grid &grid)
{
        std::ofstream openfile("/home/anton/dev/hbs/simdata/"+fname, std::ios::trunc);
        Int Nx = u.nx_, Ny = u.ny_, Nz = u.nz_;

        // Write Nsize metadata assuming Nx=Ny=Nz
        openfile << "N" << "\t" << Nx << std::endl;
        openfile << "#-----------------" << std::endl;
        
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        for (Int k=0;k<Nz;k++)
                        {
                                openfile << u(i,j,k,vi) <<  "\t" << std::endl;
                        }
                }
        }

        return 0;
                      
}
