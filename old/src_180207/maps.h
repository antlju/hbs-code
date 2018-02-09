#pragma once
#include "common.h"

void ff2sbundle(const Field4 &ff, sBundle &B, Int j, Int k)
{
        for (Int q=0;q<(Int)B.stencilSize_;q++)
        {
                Intpair jk = qtojk(q);
                for (Int i=0;i<(Int)B.nx_;i++)
                {
                        B(i,q,0) = ff(i,j+jk.first,k+jk.second,0);
                }
                //B.printpencil(q,vi);
                
        }
} // end ff2sbundle()

void ff2vbundle(const Field4 &ff, vBundle &B, Int j, Int k, Int extvi=0)
{
        if (extvi==0)
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);
                        for (Int vi=0;vi<(Int)B.nvars_;vi++)
                        {
                                for (Int i=0;i<(Int)B.nx_;i++)
                                {
                                        B(i,q,vi) = ff(i,j+jk.first,k+jk.second,vi);
                        
                                }
                                //B.printpencil(q,vi);
                        }
                }
        }
        else
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);
                        for (Int i=0;i<(Int)B.nx_;i++)
                        {
                                B(i,q,extvi) = ff(i,j+jk.first,k+jk.second,extvi);
                        
                        }
                        
                }
        }
} // end ff2vbundle()
