//This is a test of the index map in the field class
//i.e data_[vi*(zi+mz_*(yi+my_*xi))]
//it seems to work.

#include <iostream>
//#include "field.h"
#include "common.h"
#define NX 3
#define NY 3
#define NZ 3
#define NVAR 3

using namespace std;

inline size_t vfidx(Int xi, Int yi, Int zi, Int vi = 0)
{
	return vi * (NZ + 2 * NGHOST) * (NY + 2 * NGHOST) * (NX + 2 * NGHOST) +
	    (zi * (NY + 2 * NGHOST) + yi) * (NX + 2 * NGHOST) + xi;
}

int main()
{
        //Each direction has 2*NGHOST extra points
        Real u1d[(NX+2*NGHOST)*(NY+2*NGHOST)*(NZ+2*NGHOST)*NVAR] = {0.0};
        Real u4d[NVAR][(NX+2*NGHOST)][(NY+2*NGHOST)][(NZ+2*NGHOST)] = {0.0};

        Int size1d = (NX+2*NGHOST)*(NY+2*NGHOST)*(NZ+2*NGHOST)*NVAR;
        Int nvar = NVAR,mx=(NX+2*NGHOST),my=(NY+2*NGHOST),mz=(NZ+2*NGHOST);
        
        for (Int vi=0;vi<nvar;vi++){
                
                for (Int zi =NGHOST;zi<NZ+NGHOST;zi++)
                {
                        for (Int yi =NGHOST;yi<NY+NGHOST;yi++)
                        {
                                for (Int xi = NGHOST;xi<NX+NGHOST;xi++)
                                {
                                        u1d[vfidx(xi,yi,zi,vi)] = (Real)(vi+1.0);
                                        u4d[vi][xi][yi][zi] = (Real)(vi+1.0);
                                        cout << u1d[vfidx(xi,yi,zi,0)]-u4d[0][xi][yi][zi] << " ";
                                        //cout << u1d[vfidx(xi,yi,zi,vi)] << "-" << u4d[vi][xi][yi][zi] << " ";
                                }
                                cout << "x-loop" << endl;
                        }
                        cout << "y-loop" << endl;
                }
                cout << "z-loop" << endl;
                cout << "----------"<< "var: " << vi+1 << endl;
        }
        
        return 0;
}
