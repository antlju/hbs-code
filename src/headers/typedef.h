#pragma once
#include "pmesh.h"
#include "fmesh.h"
#include "fftwmesh.h"
#include "inputparams.h"

typedef fMesh<Real, NGHOSTS> Mesh;

typedef pMesh<Real, 0> Pencil;
typedef pMesh<Real, NGHOSTS> Bundle;

typedef FFTWMesh<Real> fftwMesh;
