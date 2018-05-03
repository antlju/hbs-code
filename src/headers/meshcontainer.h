#pragma once

#include "typedef.h"

/// Container class for various fields.
/// The type of constructor used will correspond to some
/// class of solver.

class MeshContainer {
public:

        /// Velocity field and applied RHS
        Mesh& u;
        Mesh& du;

        /// Runge-Kutta 4 step fields
        Mesh& k1;
        Mesh& k2;
        Mesh& k3;
        Mesh& k4;

        /// Pressure scalar field
        Mesh &p;

        /// Constructor for Kolmogorov flow Navier stokes
        
        MeshContainer(Mesh& velocity, Mesh& dvelocity, Mesh& pp) :
                u(velocity), du(dvelocity),
                k1(velocity),k2(velocity),k3(velocity),k4(velocity),
                p(pp)
        {
                
        }
        
        
}; // End class MeshContainer.
