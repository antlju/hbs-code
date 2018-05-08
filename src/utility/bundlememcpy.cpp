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
} //End pencil2ff()

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
        
} // end ff2bundle()

void ff2pencil(const Mesh &ff, Pencil &P, Int j, Int k, Int ffvi)
{
        if (P.isScalar_ == 1)
        {
                for (Int i=0;i<(Int)P.nx_;i++)
                {
                        P(i,0,0) = ff(i,j,k,ffvi);
                }
        }
        else
        {
                for (size_t vi=0;vi<P.nvars_;vi++)
                {
                        for (Int i=0;i<(Int)P.nx_;i++)
                        {
                                P(i,0,vi) = ff(i,j,k,vi);
                        }
                }
        }
} // end ff2pencil()

/// Function for copying from a bundle to a pencil.
void bundle2pencil(const Bundle &B, Pencil &P)
{
        assert((B.nx_ == P.nx_) && (B.nvars_ == P.nvars_) );
        for (size_t i=0;i<B.nx_;i++)
        {
                for (size_t vi=0;vi<B.nvars_;vi++)
                {
                        P(i,0,vi) = B(i,0,vi);
                }
        }
} //End bundle2pencil
