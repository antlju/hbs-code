#include "includes.h"

Int main()
{
        /// Time settings.
        Real dt = 0.01;
        Int maxtsteps = 100;

        /// Create parameter object and initialise parameters.
        SolverParams params;
        params.dt = dt;
        params.maxTimesteps = maxtsteps;

        /// Set grid sizes
        Real L0 = -M_PI, L1 = M_PI;
        const Int Size = 8;
        
        /// Initialise uniform 3D finite difference grid object.
        Grid<Size,Size,Size> grid(L0,L1);
        
        //ff.print();
}
