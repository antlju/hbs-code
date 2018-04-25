#include "includes.h"

/// These functions copy between the "finite-difference" main arrays (with ghost points) and
/// the FFTW3 arrays (with its own memory alignment as defined by fftw_malloc()).
/// Copy should only be done to and from the "space domain".

void mesh2fftw(const Mesh &input, fftwMesh &out, const Int vi=0)
{
        // Make sure that dimensions match.
        assert(input.nx_ == out.nx_);
        assert(input.ny_ == out.ny_);
        assert(input.nz_ == out.nz_);

        //Copy loop
        for (size_t i=0;i<input.nx_;i++)
        {
                for (size_t j=0;j<input.ny_;j++)
                {
                        for (size_t k=0;k<input.nz_;k++)
                        {
                                out.spacedom(i,j,k) = input(i,j,k,vi);
                        }
                }
        }
                
}

void fftw2mesh(const fftwMesh &input, Mesh &out, const Int vi=0)
{
        // Make sure that dimensions match.
        assert(input.nx_ == out.nx_);
        assert(input.ny_ == out.ny_);
        assert(input.nz_ == out.nz_);

        //Copy loop
        for (size_t i=0;i<input.nx_;i++)
        {
                for (size_t j=0;j<input.ny_;j++)
                {
                        for (size_t k=0;k<input.nz_;k++)
                        {
                                out(i,j,k,vi) = input.spacedom(i,j,k);
                        }
                }
        }
                
}
