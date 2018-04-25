#pragma once

#include <fftw3.h>
#include <cassert>
#include "vartypedef.h"

/// A template class of type T for storing and accessing memory
/// allocated by fftw_malloc, to be used by the FFTW3 Fourier Transform package.
/// Specifically for in-place real-to-complex 3D transforms to solve Poisson's equation.

template<class T>
class FFTWMesh {
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
        FFTWMesh(size_t Nx, size_t Ny, size_t Nz) :
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

        const fftw_complex *outPtr() const
        {
                return mem_;
        }

        T *inPtr()
        {
                return mem_[0];
        }

        const T *inPtr() const
        {
                return mem_[0];
        }

        T& spacedom(Int i, Int j, Int k)
        {

                return inPtr()[iindx(i,j,k)];
        }

        const T& spacedom(Int i, Int j, Int k) const
        {

                return inPtr()[iindx(i,j,k)];
        }

        T& freqdomRe(Int i, Int j, Int k)
        {
                return outPtr()[oindx(i,j,k)][0];
        }

        T& freqdomIm(Int i, Int j, Int k)
        {
                return outPtr()[oindx(i,j,k)][1];
        }

        const T& freqdomRe(Int i, Int j, Int k) const
        {
                return outPtr()[oindx(i,j,k)][0];
        }

        const T& freqdomIm(Int i, Int j, Int k) const
        {
                return outPtr()[oindx(i,j,k)][1];
        }


        void fill(T val)
        {
                for (size_t i=0;i<nx_;i++)
                {
                        for (size_t j=0;j<ny_;j++)
                        {
                                for (size_t k=0;k<nzh_;k++)
                                {
                                        freqdomRe(i,j,k) = val;
                                        freqdomIm(i,j,k) = val;
                                }
                        }
                }
        }

        void printspacedom()
        {
                std::cout << "========================" << std::endl;
                std::cout << "FFTWMesh 'input' print. " << std::endl;
                std::cout << "========================" << std::endl;
                for (size_t i=0;i<nx_;i++)
                {
                        for (size_t j=0;j<ny_;j++)
                        {
                                for (size_t k=0;k<nz_;k++)
                                {
                                        std::cout << spacedom(i,j,k) << "\t";
                                }
                                std::cout << std::endl;
                        }
                        std::cout << "----------------------------" << std::endl;
                }
                std::cout << std::endl;

        }

        
        /// Default destructor
        ~FFTWMesh()
        {
                fftw_free(mem_);
        }
        
};
