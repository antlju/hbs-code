#include <iostream>
#include "field.h"
#include "common.h"

using namespace std;

typedef vField<Real> RealF;

Int main()
{
        Int Nx=2,Ny=2,Nz=2;
        Int Mx=Nx+2*NGHOST,My=Ny+2*NGHOST,Mz=Nz+2*NGHOST,nvar=3;

        Real *inD = new Real[Mx*My*Mz];
        RealF uTest(nvar);
        RealF ux = uTest.subfield(1,1);
        RealF uy = uTest.subfield(2,1);
        RealF uz = uTest.subfield(3,1);
        
        return 0;
}
