#pragma once

#include "typedef.h"

/// Container class for various fields.
/// The type of constructor used will correspond to some
/// class of solver.

class MeshContainer {
public:

        /// Velocity field
        Mesh& u;
        Mesh& ustar;

        ///Runge-Kutta 3 "right hand sides" (Kim & Moin 1985)
        Mesh& RHSk; //For some k in {1,2,3}
        Mesh& RHSk_1; //For k-1
        
        /// Pressure scalar field
        Mesh &p;
        Mesh &psi;
        Mesh &gradpsi;
        fftwMesh &fftwPsi;
        
        /// Constructor for Kolmogorov flow Navier stokes
        MeshContainer(Mesh& velocity, Mesh &velstar, Mesh& rhsk, Mesh& rhsk_1, Mesh& pp, Mesh& psi_, Mesh &gradPsi, fftwMesh& fftwpsi) :
                u(velocity), ustar(velstar),
                RHSk(rhsk), RHSk_1(rhsk_1),
                p(pp), psi(psi_), gradpsi(gradPsi), fftwPsi(fftwpsi)
        {
                
        }
        
        
}; // End class MeshContainer.
