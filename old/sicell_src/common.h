#ifndef _COMMON_H_ //These things are "include guards".
#define _COMMON_H_

//some often used includes.
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cmath>

//Real space dimensions.

#define NX 8
#define NY 8
#define NZ 8

#define NGHOST 2 //Comes from the accuracy order of fin-diff scheme, in this case 4th order accuracy.

using namespace std;

//redefining things for precision considerations.
typedef double Real;
typedef int Int;

#include "field.h"
typedef vField<Real> Rfield;

#endif

