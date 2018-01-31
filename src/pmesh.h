#pragma once

#include <vector>
#include <cassert>
#include <stdexcept>
#include "vartypedefs.h" /// typedefs for double --> Real, and int --> Int.

/*
 *
 * A template "data mesh" class to be used for "pencils" (H=0) or "bundles" (H>0).
 * Based on Joonas Nättilä's cpp-toolbox: https://github.com/natj/cpp-toolbox
 *
 */

template <class T,Int Nvars, Int H=0>
class pMesh
{
public:
        /// Internal storage
        std::vector<T> mem_;

        /// Size of principal dimension.
        size_t nx_;

        /// Number of variables (i. e is it scalar pencil or vector pencil).
        size_t nvars_ = Nvars;

        size_t stencilSize_ = 4*H+1;


        /// ------ Methods --------
        /// Default constructor
        pMesh(size_t Nx) :
                nx_(Nx)
        {
                mem_.resize(Nvars*(Nx+2*H)*(4*H+1)); /// Allocate memory.
                std::fill(mem_.begin(),mem_.end(),T()); /// Fill with zeros.
        }

        /// Fill with some constant
        void fillUp(const T& val)
        {
                std::fill(mem_.begin(),mem_.end(),val);
        }
        
        /// Clear internal storage (overwrite with zeros)
        void clear()
        {
                std::fill(mem_.begin(),mem_.end(), T() ); 
        }
        
        /// Validate compatible sizes of lhs and rhs.
        void validateDims(const pMesh &rhs)
        {
                if(this->nx_ != rhs.nx_)
                        throw std::range_error ("lhs and rhs dimensions do no match!");
        }

        /// Indexing for pencils and bundles
        size_t pencil_index(Int i,Int vi=0) const
        {
                Int pindx = vi*nx_+i;
                return pindx;
        }

        size_t bundle_index(Int i, Int q, Int vi=0) const
        {
                
                Int bindx = vi*(nx_+2*H)+(i+H);
                return bindx;
        }
        
        /// ------ Overloads --------
        /// Element accessing (This will contain proper maps later!)
        
        ///           pencil element access:
        T& operator()(Int i,Int vi=1)
        {
                return mem_[pencil_index(i,vi)];
        }
        const T& operator()(Int i,Int vi=1) const
        {
                return mem_[pencil_index(i,vi)];
        }

        ///           bundle element access:
        T& operator()(Int i)
        {
                return mem_[i];
        }
        const T& operator()(Int i) const
        {
                return mem_[i];
        }
        
        /// Arithmetics.
        pMesh& operator=(const pMesh& rhs);

        pMesh& operator=(const T& rhs);

        pMesh& operator+=(const pMesh& rhs);

        pMesh& operator-=(const pMesh& rhs);

        pMesh& operator*=(const T& rhs);

        pMesh& operator/=(const T& rhs);


}; /// End class pMesh

/// -------------------
///      Arithmetics operator overloading
/// -------------------
template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator=(const pMesh<T,Nvars,H>& rhs)
{
        validateDims(rhs);

        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] = rhs.mem_[i];
        }

        return *this;
}

template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator=(const T& rhs)
{
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] = rhs;
        }
}

template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator+=(const pMesh<T,Nvars,H>& rhs)
{
        validateDims(rhs);
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] += rhs.mem_[i];
        }
        
}

template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator-=(const pMesh<T,Nvars,H>& rhs)
{
        validateDims(rhs);
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] -= rhs.mem_[i];
        }
        return *this;  
}

template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator*=(const T& rhs)
{
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] *= rhs;
        }
        return *this;
}

template <class T, Int Nvars, Int H>
pMesh<T,Nvars,H>& pMesh<T,Nvars,H>::operator/=(const T& rhs)
{
        for (size_t i=0;i<this->mem_.size();i++)
        {
                this->mem_[i] /= rhs;
        }
        return *this;
}

///-------------------------------------------------- 
///    Array arithmetics, ie pencilA+pencilB = pencilC etc.
///--------------------------------------------------
template <class T, Int Nvars, Int H>
inline pMesh<T,Nvars,H> operator+(pMesh<T,Nvars,H> lhs, const pMesh<T,Nvars,H>& rhs)
{
        lhs += rhs;
        return lhs;
}

template <class T, Int Nvars, Int H>
inline pMesh<T,Nvars,H> operator-(pMesh<T,Nvars,H> lhs, const pMesh<T,Nvars,H>& rhs)
{
        lhs -= rhs;
        return lhs;
}

///-------------------------------------------------- 
///    Single value operators
///--------------------------------------------------
template <class T, Int Nvars, Int H>
inline pMesh<T,Nvars,H> operator*(pMesh<T,Nvars,H> lhs, const T& rhs)
{
        lhs *= rhs;
        return lhs;
}

template <class T, Int Nvars, Int H>
inline pMesh<T,Nvars,H> operator/(pMesh<T,Nvars,H> lhs, const T& rhs)
{
        lhs /= rhs;
        return lhs;
}
