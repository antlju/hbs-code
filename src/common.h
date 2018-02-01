#pragma once

#include <iostream>
#include <cmath>

#define NGHOSTS 2
#define NX 4
#define NY 4
#define NZ 4

#include "vartypedefs.h"
#include "pmesh.h"
#include "field.h"

typedef Field<Real,NGHOSTS> Field4;

typedef pMesh<Real,1,0> sPencil;
typedef pMesh<Real,3,0> vPencil;

typedef pMesh<Real,1,NGHOSTS> sBundle;
typedef pMesh<Real,3,NGHOSTS> vBundle;

#include "diffops.h"
