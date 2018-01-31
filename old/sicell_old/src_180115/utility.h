#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "common.h"

inline size_t vfidx(Int xi, Int yi, Int zi, Int vi = 0)
{
	return vi * (NZ + 2 * NGHOST) * (NY + 2 * NGHOST) * (NX + 2 * NGHOST) +
	    (zi * (NY + 2 * NGHOST) + yi) * (NZ + 2 * NGHOST) + xi;
}

inline size_t vfidx1D(Int i)
{
	return i+NGHOST;
}


Real *linspace(Int size, Real start, Real end)
{
        Real dx = (end-start)/size;
        Real *axis = new Real[size];
        for (int i=0;i<size;i++)
        {
                axis[i] = start+dx*i;
        }
        return axis;
}

#endif
