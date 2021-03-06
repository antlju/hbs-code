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

class PBContainer {
public:
	Bundle &uBndl; // velocity, velocity* and pressure and psi bundles
	Bundle &ustarBndl; 
	Bundle &pBndl;
	Bundle &psiBndl;

	Pencil &uPncl; // old u as pencil
	Pencil &ustarPncl;
	Pencil &dvPncl;// vector result of diffop
	Pencil &dsPncl; //scalar result of diffop
	Pencil &fPncl; //Force
	Pencil &rhskPncl; // RHSk RK substep pencils
	Pencil &rhsk_1Pncl;
	
PBContainer(Bundle &u, Bundle &ustar, Bundle &p, Bundle &psi,
	    Pencil &up, Pencil &ustarp, Pencil &dv, Pencil &ds,
	    Pencil &f, Pencil &rhsk, Pencil &rhsk_1) :
	uBndl(u), ustarBndl(ustar), pBndl(p), psiBndl(psi),
		uPncl(up), ustarPncl(ustarp), dvPncl(dv), dsPncl(ds),
		fPncl(f), rhskPncl(rhsk), rhsk_1Pncl(rhsk_1)
		
	{
		
	}
		
}; // End class PBContainer
