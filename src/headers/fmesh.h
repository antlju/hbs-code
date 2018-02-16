#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include "vartypedef.h"

//A class of type T with NG number of ghost points.
//Based on Joonas Nättilä's cpp-toolbox: https://github.com/natj/cpp-toolbox
//For example NG=2 for 4th order accuracy finite difference.

template<class T, Int NG=0> 
class fMesh {
public:
        /// Internal storage
        std::vector<T> mem_;

        /// grid size along x (excluding ghosts)
        size_t nx_;

        /// grid size along y (excluding ghosts)
        size_t ny_;

        /// grid size along z (excluding ghosts)
        size_t nz_;

        /// number of variables (i.e if scalar field: nvar_=1, or a vector field nvar_=3 etc.)
        size_t nvar_;

        Int isScalar_ = 0;
        
        size_t ng_ = NG;
        
        /// Internal indexing function taking ghostpoints into account. I.e 4D->1D map
        size_t indx(Int i, Int j, Int k, Int vi) const
        {
                Int indx = vi*(nz_+2*NG)*(ny_+2*NG)*(nx_+2*NG)
                        +(i+NG)+(ny_+2*NG)*((j+NG)+(nz_+2*NG)*(k+NG)); 
                return indx;
        }

        /// This index map is used with operator overloading
        T& operator()(Int i, Int j, Int k, Int vi=0)
        {
                return mem_[ indx(i,j,k,vi) ];
        }

        const T& operator()(Int i, Int j, Int k, Int vi=0) const
        {
                return mem_[ indx(i,j,k,vi) ];
        }
        

        /// Default constructor
        fMesh(size_t Nx, size_t Ny, size_t Nz, size_t Nvars) :
                nx_(Nx), ny_(Ny), nz_(Nz), nvar_(Nvars)
        {
                mem_.resize( Nvars*(Nx+2*NG)*(Ny+2*NG)*(Nz+2*NG) ); //Allocate memory
                std::fill(mem_.begin(),mem_.end(), T() ); //Fill with zeros of type T.
                if (Nvars == 1)
                        isScalar_ = 1;
        }

        /// Validate compatible sizes of lhs and rhs.
        void validateDims(const fMesh &rhs)
        {
                if( (this->nx_ != rhs.nx_)
                    || (this->ny_ != rhs.ny_)
                    || (this->nz_ != rhs.nz_) )
                        throw std::range_error ("lhs and rhs dimensions do no match!");
        }
        
        /// Arithmetics.
        fMesh& operator*=(const T& rhs);
        fMesh& operator=(const fMesh& rhs);
        fMesh& operator+=(const fMesh& rhs);
        
                      
        
}; //End class

template <class T, Int NG>
fMesh<T,NG>& fMesh<T,NG>::operator=(const fMesh<T,NG>& rhs)
{
        validateDims(rhs);

        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] = rhs.mem_[i];
        }

        return *this;
}

template <class T, Int NG>
fMesh<T,NG>& fMesh<T,NG>::operator*=(const T& rhs)
{
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] *= rhs;
        }
        return *this;
}

template <class T, Int NG>
fMesh<T,NG>& fMesh<T,NG>::operator+=(const fMesh<T,NG>& rhs)
{
        validateDims(rhs);
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] += rhs.mem_[i];
        }

        return *this;
}

///-------------------------------------------------- 
///    Array arithmetics, ie meshA+meshB = meshC etc.
///--------------------------------

template <class T, Int NG>
inline fMesh<T,NG>& operator*(fMesh<T,NG> lhs, const T& rhs)
{
        lhs *= rhs;
        return lhs;
}

template <class T, Int H>
inline fMesh<T, H> operator+(fMesh<T, H> lhs, const fMesh<T, H>& rhs)
{
        lhs += rhs;
        return lhs;
}

