#pragma once
//#include "common.h"
//#include "qjkmap.h"

void ff2sbundle(const Field4 &ff, sBundle &B, Int j, Int k);
void spencil2ff(const sPencil &P, Field4 &ff, Int j, Int k,Int vi=0);
void ff2vbundle(const Field4 &ff, vBundle &B, Int j, Int k, Int extvi=0);
