#pragma once

#include <fftw3.h>
#include <cassert>
#include "vartypedef.h"

/// A template class of type T for storing and accessing memory
/// allocated by fftw_malloc, to be used by the FFTW3 Fourier Transform package.
/// Specifically for use in in-place real-to-complex 3D transforms for solving
/// Poisson's equation by Fourier transform.

template<class T>
class fftwMesh {
public:
        /// Internal main storage
        fftw_complex *mem_;

        /// input pointer
        T *input_ = mem_[0];

        /// output (that is Fourier transformed data)
        fftw_complex *output_ = mem_;

        /// grid size along x
        size_t nx_;

        /// grid size along y
        size_t ny_;

        /// grid size along z
        size_t nz_;

        /// Slightly more than half z-size (Nzh = Nz/2+1)
        size_t nzh_;
        
        /// Default constructor
        fftwMesh(size_t Nx, size_t Ny, size_t Nz) :
                nx_(Nx), ny_(Ny), nz_(Nz), nzh_(Nz/2+1)
        {
                mem_ = (fftw_complex*) fftw_malloc(
                        sizeof(fftw_complex) * Nx*Ny*(Nz/2+1));
        }

        /// Internal indexing for "input" (real data) array
        size_t iindx(size_t i, size_t j, size_t k) const
        {
                assert(i<nx_ && j < ny_ && k < nz_);
                return k+(nz_+2)*(j+ny_*i);
        }

        /// Internal indexing for "input" (real data) array
        size_t oindx(size_t i, size_t j, size_t k) const
        {
                assert(i<nx_ && j < ny_ && k < nzh_);
                return k+(nzh_)*(j+ny_*i);
        }

        fftw_complex *outPtr()
        {
                return mem_;
        }

        T *inPtr()
        {
                return mem_[0];
        }

        T& in(Int i, Int j, Int k)
        {

                return inPtr()[iindx(i,j,k)];
        }

        T& outRe(Int i, Int j, Int k)
        {
                return outPtr()[oindx(i,j,k)][0];
        }

        T& outIm(Int i, Int j, Int k)
        {
                return outPtr()[oindx(i,j,k)][1];
        }

        
        /// Default destructor
        ~fftwMesh()
        {
                fftw_free(mem_);
        }
        
};
