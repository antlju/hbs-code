#include "common.h"
#include "diffops.h"
#include "maps.h"

/// "OLD STYLE" CURL THIS WORKS!

void bundlecurl(vBundle &B, sPencil &xp, sPencil &yp, sPencil &zp, const Real xfac, const Real yfac,const Real zfac)
{
        //B.printpencil(0,0);
        Real d1,d2;
        for (size_t i=0;i<B.nx_;i++)
        {
               
                //Compute omega_1 = d_2u_3-d_3u_2
                d1 = yfac * der1(B(i,7,2),B(i,3,2), //refer to qjkmap.h for q->j vals.
                                     B(i,1,2),B(i,5,2));
                
                d2 = zfac * der1(B(i,6,1),B(i,2,1), //refer to qjkmap.h for q->k vals.
                                     B(i,4,1),B(i,8,1));

                //std::cout << d1 << "\t" << d2 << "\t" << d1-d2 << std::endl;
                xp(i) = d1-d2;

                //Compute omega_2 = d_3u_1-d_1u_3
                d1 = zfac * der1(B(i,6,0),B(i,2,0), //refer to qjkmap.h for q->k vals.
                                     B(i,4,0),B(i,8,0));

                d2 = xfac * der1(B(i-2,0,2),B(i-1,0,2),
                                     B(i+1,0,2),B(i+2,0,2));

                yp(i) = d1-d2;

                //Compute omega_3 = d_1u_2-d_2u_1
                d1 = xfac * der1(B(i-2,0,1),B(i-1,0,1),
                                     B(i+1,0,1),B(i+2,0,1));
                
                d2 = yfac * der1(B(i,7,0),B(i,3,0), //refer to qjkmap.h for q->j vals.
                                     B(i,1,0),B(i,5,0));

                zp(i) = d1-d2;
                //std::cout << xp(i) << "\t" << yp(i) << "\t" << zp(i) << std::endl;
        }
}


//cu = curl(u), it's a map curl:vector field -> vector field.
void curl(const Field4 &u, Field4 &cu, const Real xfactor, const Real yfactor, const Real zfactor)
{
        vBundle Dbundle((Int)u.nx_);
        sPencil xpencil((Int)u.nx_);
        sPencil ypencil((Int)u.nx_);
        sPencil zpencil((Int)u.nx_);
        
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {
                        /// Copies from big mem properly to all bundle components!
                        /// This I have checked!
                        ff2vbundle(u,Dbundle,j,k,0);
                        ff2vbundle(u,Dbundle,j,k,1);
                        ff2vbundle(u,Dbundle,j,k,2);

                        //Dbundle.printpencil(0,0);
                        //Dbundle.printpencil(0,1);
                        //Dbundle.printpencil(0,2);
                        /// Compute the curl on a bundle
                        bundlecurl(Dbundle,
                                   xpencil,ypencil,zpencil,
                                   xfactor,yfactor,zfactor); 

                        
                        /// Copy resulting pencils to big mem
                        spencil2ff(xpencil,cu,j,k,0);
                        spencil2ff(ypencil,cu,j,k,1);
                        spencil2ff(zpencil,cu,j,k,2);
                }
        }
}
