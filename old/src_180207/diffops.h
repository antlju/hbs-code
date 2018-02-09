#pragma once
#include "common.h"

Real der1(const Real m2h, const Real m1h, const Real p1h, const Real p2h);
Real der2(const Real m2h, const Real m1h, const Real mid, const Real p1h, const Real p2h);

void grad(const sBundle &f, vPencil &g);
