#ifndef _FIELD_H_
#define _FIELD_H_

#include "common.h"

// An N-dimensional field is
// written u(x) = (u_1(x),u_2(x),...,u_N(x)),
// where x = x_1,x_2,...,x_N. Here the case is N=3 and cartesian labels x,y,z.
// Data is allocated contiguous as a long one-dimensional array. 
// Thus we need to allocate Mx*My*Mz*3 of space for a function u(x) = (u_x(x),u_y(x),u_z(x)).


//I sort of understand the stride. It's the separation of the axes in memory.
inline size_t ystride() //y-elements are separated by a line of x-elements.
{
        return (NX+2*NGHOST);
}

inline size_t zstride()//z-elements are separated by xy-plane.
{
        return (NY+2*NGHOST) *ystride();
}

//But the point is to use xystride() multiplied by variable index to map to the desired location
//in allocated memory. That is xyzstride() returns MX*MY*MZ, where Mi = Ni+2*NGHOST.
inline size_t xyzstride() //Here we stride an entire space function (i.e one "variable")
{
        return (NZ+2*NGHOST) *zstride();
}

//A function for calculating memory needed for a spacefield of nvar variables.
inline size_t fmemsize(Int nvar)
{
        return nvar * xyzstride() * sizeof(Real);
}

// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

//Main tensor field class (can be scalar, vector, matrix etc)

template<class T>
class vField {
	T *m_mem;
	Int m_nvar;
        vField(T *mem, Int nvar);
public:
        //Constructors
	vField();
	vField(Int nvar);

        //Operator overloading
        inline T & operator[](const Int i);
        inline const T & operator[](const Int i) const;
        
        //Subfield construction (pointers to variables, eg u_x(x,y,z) == u.subfield(1,1)
	vField<T> subfield(Int vi, Int nvar); //we index components by vi = 0,1,2,...,m_nvar-1.

        Int nvar(){return m_nvar;}
        //Destructor 
        //~vField();
};

//------------------- Member functions ---------------

//--------- Constructors
template<class T>
vField<T>::vField(T *mem, Int nvar) : m_mem(mem), m_nvar(nvar) { }

template<class T>
vField<T>::vField(): m_mem(NULL), m_nvar(0) { }

template<class T>
vField<T>::vField(Int nvar): m_nvar(nvar)
{
        Int arraysize = xyzstride()*nvar;
        
        m_mem = new T[xyzstride()*nvar];
        //m_mem[xyzstride()*1] = 1;
        //m_mem[xyzstride()*2] = 2;
        for (Int i=0;i<arraysize;i++)
        {
                m_mem[i] = 0;
        }
}

//--------- Operator overloading
template<class T>
inline T & vField<T>::operator[](const Int i)
{
        return m_mem[i];
}

//--------- Operator overloading
template<class T>
inline const T & vField<T>::operator[](const Int i) const
{
        return m_mem[i];
}

//--------- Subfield construction (pointers to variables, eg u_x(x,y,z) == u.subfield(1,1)
template<class T>
vField<T> vField<T>::subfield(Int vi, Int nvar)
{
        assert(vi >= 0 && vi<=m_nvar-1); //we index components by vi = 0,1,2,...,m_nvar-1.
        return vField(m_mem + xyzstride() * vi, nvar);
}

/*
//--------- Destructor
template<class T>
vField<T>::~vField()
{
        if (m_mem != NULL) delete[] m_mem;
}
*/           


// ----- functions for getting indices, i.e 4D <-> 1D map.
inline size_t vfidx(Int xi, Int yi, Int zi, Int vi = 0)
{
	return vi * (NZ + 2 * NGHOST) * (NY + 2 * NGHOST) * (NX + 2 * NGHOST) +
	    (zi * (NY + 2 * NGHOST) + yi) * (NZ + 2 * NGHOST) + xi;
}

inline size_t vfidx1D(Int i)
{
	return i+NGHOST;
}

// ----- functions for getting indices, i.e 4D <-> 1D map.
inline size_t mat2lin(Int xi, Int yi)
{
	return (3*yi)+xi;
}



        
#endif //end #define _FIELD_H_

