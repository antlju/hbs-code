//This is a test of the index map in the field class
//i.e data_[vi*(zi+mz_*(yi+my_*xi))]
//it seems to work.

#include <iostream>
#include "field.h"
#include "common.h"

using namespace std;

typedef Field<Real> Rfield;

Int main()
{
        Int Nx=3,Ny=3,Nz=3;
        Int Mx=Nx+NGHOST,My=Ny+NGHOST,Mz=Nz+NGHOST,nvar=3,nvar3d=1;
        Rfield maptest(Nx,Ny,Nz,nvar3d); //should be 3x3x3x1 tensor.
        Real maptArr[3*3*3*1];
        for (Int iterator = 0;iterator<maptest.totSize();iterator++)
                maptArr[iterator] = static_cast<double>(iterator);

        Int iter = 0;
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        for (Int k=0;k<Nz;k++)
                        {
                                maptest.set(i,j,k,1,static_cast<double>(iter));
                                iter++;
                        }
                }
        }

        iter = 0;
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        for (Int k=0;k<Nz;k++)
                        {
                                cout << maptArr[iter]-maptest.get(i,j,k,1) << endl;
                                iter++;
                        }
                }
        }


        

        
        
        
        
        return 0;
}
