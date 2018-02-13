#pragma once
void linspace(Pencil &u, Real start, Real end, Real dx);
void apply_pbc(Mesh &u);
void pencil2ff(const Pencil &P, Mesh &ff, Int j, Int k, Int ffvi=0);
void ff2bundle(const Mesh &ff, Bundle &B, Int j, Int k, Int ffvi=0);
void printfield(const Mesh &u);
