#include "includes.h"

Int writeStatsToFile(const std::string fname, const MeshContainer &meshCntr, const Stats &stats, const Grid &grid, const Int stepNo)
{
	std::ofstream openfile("../simdata/"+fname, std::ios::app); //Append option!

	Real volsize = 3*meshCntr.u.nx_*meshCntr.u.ny_*meshCntr.u.nz_;
	
	Real avgomega2 = stats.omega2/volsize;
	//Real avgP = stats.P/volsize;
	//Real avgP2 = stats.P2/volsize;
	Real avgE = stats.energy/volsize;
	//std::cout << "stats.urms_step" << stats.urms_step << std::endl;
	//Real urms = sqrt(stats.urms_step/volsize);
	Real urms_mesh = sqrt(stats.urms_mesh/volsize);

	//Real Re = urms_mesh*2*M_PI/(10.0);
	//openfile << stepNo << ":\t" << avgomega2 << "\t" << avgP << "\t" << avgP2 << "\t" << stats.umax << "\t" << avgE << std::endl;
	//openfile << stepNo << std::setw(8) << std::left << avgomega2 << std::setw(8) << std::left << urms << std::setw(8) << std::left << urms_mesh << std::setw(8) << std::left << stats.umax << std::setw(8) << std::left << avgE << std::endl;
		openfile << stepNo << "    " << avgomega2 << "    " << urms_mesh  << "    " << stats.umax << "    " << avgE << std::endl;
	
	return 0;
}

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
        std::ofstream openfile("../simdata/"+fname, std::ios::trunc);
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
        std::ofstream openfile("../simdata/"+fname, std::ios::trunc);
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

	std::cout << "wrote " << fname << std::endl;
        return 0;
                      
}



/* // This old one writes to a new file each time 
Int writeStatsToFile(const std::string fname, const MeshContainer &meshCntr, const Stats &stats, const Grid &grid)
{
	std::ofstream openfile("../simdata/"+fname, std::ios::trunc);

	Real volsize = meshCntr.u.nx_*meshCntr.u.ny_*meshCntr.u.nz_;
	
	Real avgomega2 = stats.omega2/volsize;
	Real avgP = stats.P/volsize;
	Real avgP2 = stats.P2/volsize;
	
	openfile << "<Omega^2> \t" << "<P> \t" << "<P^2> \t" << "umax" << std::endl;
	openfile << avgomega2 << "\t" << avgP << "\t" << avgP2 << "\t" << stats.umax << std::endl;

	return 0;
}
*/
