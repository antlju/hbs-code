// An N-dimensional field is
// written u(x) = (u_1(x),u_2(x),...,u_N(x)),
// where x = x_1,x_2,...,x_N. Here the case is N=3 and cartesian labels x,y,z.
// Data is allocated contiguous as a long one-dimensional array. 
// Thus we need to allocate Mx*My*Mz*3 of space for a function u(x) = (u_x(x),u_y(x),u_z(x)).

#include <assert.h>
#include <iostream>
#include "common.h"

//I don't understand why the strides are defined this way yet.
inline size_t ystride()
{
        return (NX+2*NGHOST);
}

inline size_t zstride()
{
        return (NY+2*NGHOST) *ystride();
}

//But the point is to use xystride() multiplied by variable index to map to the desired location
//in allocated memory. That is xyzstride() returns MX*MY*MZ, where Mi = Ni+2*NGHOST.
//
inline size_t xyzstride()
{
        return (NZ+2*NGHOST) *zstride();
}

//A function for calculating memory space needed for a spacefield of nvar variables.
inline size_t fmemsize(Int nvar)
{
        return nvar * xyzstride() * sizeof(Real);
}

template <class T>
T testPtrF(T *test){
        return test[0];
}

//Main field class 
template<class T>
class vField {
	T *m_mem;
	Int m_nvar;
	vField(T *mem, Int nvar): m_mem(mem), m_nvar(nvar) { std::cout << m_mem[0] << std::endl;}
public:
	vField(): m_mem(NULL), m_nvar(0) { }
	vField(Int nvar): m_nvar(nvar) {
		m_mem = new T[xyzstride()*nvar];
                m_mem[xyzstride()*1] = 1;
                m_mem[xyzstride()*2] = 2;
	}

	vField<T> subfield(Int vi, Int nvar) {
                assert(vi > 0 && vi<=m_nvar); //we index components by vi = 1,2,3,...,m_nvar.
		return vField(m_mem + xyzstride() * (vi-1), nvar);
        }
};

        


