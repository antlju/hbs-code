#pragma once

#include <iostream>
#include <cmath>

#define NGHOSTS 2
#define NX 8
#define NY 8
#define NZ 8

#include "vartypedefs.h"
#include "field.h"

typedef Field<Real,NGHOSTS> Rfield;
