#include "includes.h"
#include "qjkmap.h"

void pencil2ff(const Pencil &P, Mesh &ff, Int j, Int k,Int ffvi)
{
        if (P.isScalar_ == 1)
        {
                for (size_t i=0;i<ff.nx_;i++)
                {
                        ff(i,j,k,ffvi) = P(i,0,0);
                }
        }
        else
        {
                for (size_t vi=0;vi<ff.nvar_;vi++)
                {
                        for (size_t i=0;i<ff.nx_;i++)
                        {
                                ff(i,j,k,vi) = P(i,0,vi);
                        }
                }
        }
}

void ff2bundle(const Mesh &ff, Bundle &B, Int j, Int k, Int ffvi)
{
        if (B.isScalar_ == 1)
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);

                        for (Int i=0;i<(Int)B.nx_+2*(Int)ff.ng_;i++)
                        {
                                B(i-(Int)ff.ng_,q,0) =
                                        ff(i-(Int)ff.ng_,j+jk.first,k+jk.second,ffvi);
                        
                        }
                       
                        
                }
        }
        else
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);
                        for (size_t vi=0;vi<B.nvars_;vi++)
                        {
                                for (Int i=0;i<(Int)B.nx_+2*(Int)ff.ng_;i++)
                                {
                                        B(i-(Int)ff.ng_,q,vi) =
                                                ff(i-(Int)ff.ng_,j+jk.first,k+jk.second,vi);
                        
                                }
                        }
                        
                }
        }
        
} // end ff2vbundle()
