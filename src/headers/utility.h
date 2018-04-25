#pragma once
void linspace(Pencil &u, Real start, Real end, Real dx);
void apply_pbc(Mesh &u);

void pencil2ff(const Pencil &P, Mesh &ff, Int j, Int k, Int ffvi=0);
void ff2bundle(const Mesh &ff, Bundle &B, Int j, Int k, Int ffvi=0);
void bundle2pencil(const Bundle &B, Pencil &P);
        
void printfield(const Mesh &u);

void mesh2fftw(const Mesh &input, fftwMesh &out, const Int vi);
void fftw2mesh(const fftwMesh &input, Mesh &out, const Int vi);
