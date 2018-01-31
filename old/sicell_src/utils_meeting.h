#ifndef _UTILS_H_
#define _UTILS_H_

void periodic_bc(Rfield U)
{
        for (Int vi=0;vi<U.nvar();vi++)
        {
                for (Int zi=NGHOST;zi<NZ+NGHOST;zi++)
                {
                        for (Int yi=NGHOST;yi<NY+NGHOST;yi++)
                        {
                                for (Int xi=NGHOST;xi<NX+NGHOST;xi++)
                                {
                                        //x PBCs
                                        
                                        U[vfidx(0,yi,zi,vi)] = U[vfidx(NX+NGHOST-3,yi,zi,vi)];
                                        U[vfidx(1,yi,zi,vi)] = U[vfidx(NX+NGHOST-2,yi,zi,vi)];
                                        U[vfidx(NX+NGHOST,yi,zi,vi)] = U[vfidx(NGHOST+1,yi,zi,vi)];
                                        U[vfidx(NX+NGHOST+1,yi,zi,vi)] = U[vfidx(NGHOST+2,yi,zi,vi)];

                                        //y PBCs
                                        U[vfidx(xi,0,zi,vi)] = U[vfidx(xi,NY+NGHOST-3,zi,vi)];
                                        U[vfidx(xi,1,zi,vi)] = U[vfidx(xi,NY+NGHOST-2,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST,zi,vi)] = U[vfidx(xi,NGHOST+1,zi,vi)];
                                        U[vfidx(xi,NY+NGHOST+1,zi,vi)] = U[vfidx(xi,NGHOST+2,zi,vi)];
                                
                                        //z PBCs
                                        U[vfidx(xi,yi,0,vi)] = U[vfidx(xi,yi,NZ+NGHOST-3,vi)];
                                        U[vfidx(xi,yi,1,vi)] = U[vfidx(xi,yi,NZ+NGHOST-2,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST,vi)] = U[vfidx(xi,yi,NGHOST+1,vi)];
                                        U[vfidx(xi,yi,NZ+NGHOST+1,vi)] = U[vfidx(xi,yi,NGHOST+2,vi)];
                              
                                }    
                        }    
                }
        }
}

Real *linspace(Int size, Real start, Real end)
{
        Real dx = (end-start)/(size-1);
        Real *axis = new Real[size];
        for (Int i=0;i<size;i++)
        {
                axis[i] = start+dx*i;
                cout << axis[i] << endl;
        }
        cout << endl;
        return axis;
}

#endif
