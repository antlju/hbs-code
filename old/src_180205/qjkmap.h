#pragma once

#include "common.h"
#include <map>
#include <utility>

typedef std::map<Int, Int> iiMap;

Intpair qtojk(Int q)
{
        iiMap qtoj;
        iiMap qtok;
        
        qtoj[0] = 0;
        qtoj[1] = 1;
        qtoj[2] = 0;
        qtoj[3] = -1;
        qtoj[4] = 0;
        qtoj[5] = 2;
        qtoj[6] = 0;
        qtoj[7] = -2;
        qtoj[8] = 0;

        qtok[0] = 0;
        qtok[1] = 0;
        qtok[2] = -1;
        qtok[3] = 0;
        qtok[4] = 1;
        qtok[5] = 0;
        qtok[6] = -2;
        qtok[7] = 0;
        qtok[8] = 2;

        auto jmap = qtoj.find(q);
        auto kmap = qtok.find(q);
        Int j = jmap->second;
        Int k = kmap->second;

        Intpair ret(j,k);
        return ret;
        
        //return ret;
        //std::cout << j << "\t" << k << std::endl;

        
        
}
