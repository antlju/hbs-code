#include "common.h"
#include "qjkmap.h"

void ff2sbundle(const Field4 &ff, sBundle &B, Int j, Int k)
{
        for (Int q=0;q<(Int)B.stencilSize_;q++)
        {
                Intpair jk = qtojk(q);
                for (Int i=0;i<(Int)B.nx_+2*(Int)ff.ng_;i++)
                {
                        //std::cout << "i: " << i << " q: " << q << std::endl;
                        B(i-(Int)ff.ng_,q,0) = ff(i-(Int)ff.ng_,j+jk.first,k+jk.second,0);
                }
                B.printpencil(q,0);
                
        }
        //B.printpencil(0,0);
} // end ff2sbundle()

void spencil2ff(const sPencil &P, Field4 &ff, Int j, Int k,Int vi=0)
{
        for (size_t i=0;i<ff.nx_;i++)
        {
                ff(i,j,k,vi) = P(i,0);
        }
}

void ff2vbundle(const Field4 &ff, vBundle &B, Int j, Int k, Int extvi)
{
        
        for (Int q=0;q<(Int)B.stencilSize_;q++)
        {
                Intpair jk = qtojk(q);
                for (Int i=0;i<(Int)B.nx_+2*(Int)ff.ng_;i++)
                {
                        B(i-(Int)ff.ng_,q,extvi) =
                                ff(i-(Int)ff.ng_,j+jk.first,k+jk.second,extvi);
                        
                }
                        
        }
        
} // end ff2vbundle()
