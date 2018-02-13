#pragma once
#include <map>

#include "pmesh.h"
#include "fmesh.h"
#include "inputparams.h"

typedef std::pair<Int, Int> Intpair;
typedef std::map<Int, Int> iiMap;

typedef fMesh<Real, NGHOSTS> Mesh;

typedef pMesh<Real, 0> Pencil;
typedef pMesh<Real, NGHOSTS> Bundle;
