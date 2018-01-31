#include "field.h"

//Constructors
template<class T>
vField<T>::vField(T *mem, Int nvar) : m_mem(mem), m_nvar(nvar) { }

template<class T>
vField<T>::vField(): m_mem(NULL), m_nvar(0) { }

template<class T>
vField<T>::vField(Int nvar): m_nvar(nvar)
{
		m_mem = new T[xyzstride()*nvar];
                m_mem[xyzstride()*1] = 1;
                m_mem[xyzstride()*2] = 2;
}

template<class T>
vField<T> vField<T>::subfield(Int vi, Int nvar)
{
        assert(vi > 0 && vi<=m_nvar); //we index components by vi = 1,2,3,...,m_nvar.
        return vField(m_mem + xyzstride() * (vi-1), nvar);
}

template<class T>
vField<T>::~vField()
{
        if (m_mem != NULL) delete[] m_mem;
}
              
