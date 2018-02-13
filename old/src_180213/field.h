#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include "vartypedefs.h"

//A class of type T with NG number of ghost points.
//Based on Joonas Nättilä's cpp-toolbox: https://github.com/natj/cpp-toolbox
//For example NG=2 for 4th order accuracy finite difference.
template<class T, Int NG=0> 
class Field {
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
        Field(size_t Nx, size_t Ny, size_t Nz, size_t vi) :
                nx_(Nx), ny_(Ny), nz_(Nz), nvar_(vi)
        {
                mem_.resize( vi*(Nx+2*NG)*(Ny+2*NG)*(Nz+2*NG) ); //Allocate memory
                std::fill(mem_.begin(),mem_.end(), T() ); //Fill with zeros of type T.
        }

        /// Arithmetics.
        Field& operator*=(const T& rhs);
        
}; //End class

template <class T, Int NG>
Field<T,NG>& Field<T,NG>::operator*=(const T& rhs)
{
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] *= rhs;
        }
        return *this;
}
/*
template <class T, Int NG>
inline Field<T,NG>& Field<T,NG>::operator*(Field<T,NG> lhs, const T& rhs)
{
        lhs *= rhs;
        return lhs;
}
*/
