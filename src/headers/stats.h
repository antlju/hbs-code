#pragma once

#include "vartypedef.h"
#include "typedef.h"

/// Statistics calculation and storage class

class Stats {
public:
        /// max(abs(u)), The maximum absolute velocity
        Real umax;
        Real umax_old; /// 

        /// Calculate vector pencil max(abs(u))
        void calc_pencil_umax(Pencil &u)
        {
                umax_old = umax;
                for (size_t i=0;i<u.nx_;i++)
                {
                        Real absu = fabs(u(i,0,0)+u(i,0,1)+u(i,0,2));
                        if (absu > umax)
                                umax = absu;
                }

            
        }
}; //End class Stats
