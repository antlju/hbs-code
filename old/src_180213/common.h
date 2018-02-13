#pragma once

#include <iostream>
#include <cmath>

#define NGHOSTS 2
#define NX 8
#define NY 8
#define NZ 8

#include "pmesh.h"
#include "field.h"
#include "vartypedefs.h"

typedef Field<Real, NGHOSTS> Field4;

typedef pMesh<Real, 1, 0> sPencil;
typedef pMesh<Real, 3, 0> vPencil;

typedef pMesh<Real, 1, NGHOSTS> sBundle;
typedef pMesh<Real, 3, NGHOSTS> vBundle;


