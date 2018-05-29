import numpy as np
from matplotlib import pyplot as plt

def get_mesh_from_file(fname,nn,kf):
    fname0 = fname  

    u0 = np.loadtxt(fname0,skiprows=1)
    N = nn 
    Nx = N
    Ny = N
    Nz = N
    
    U0 = np.zeros((N,N,N));
    
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                U0[i][j][k] = u0[(Ny*i+j)*Nz+k]

    return U0
                
        

def calc_urms(umesh):
    urms = 0
    N = len(umesh[:][0][0])
    print("N = %i",N)
    
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                urms = urms+np.power(umesh[i][j][k],2)

    urms = np.sqrt(urms/(N*N*N))

    print urms
    return urms


def calc_urms_list():

    fname0 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3400.dat"
    fname1 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3500.dat"
    fname2 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3600.dat"
    fname3 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3700.dat"
    fname4 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3800.dat"
    fname5 = "/home/anton/code/hbs-code/simdata/steadyState_test_kf_1_component_0_N_256_stepNo_3900.dat"
    
    fnlist = [fname0,fname1,fname2,fname3,fname4,fname5]

    urms_list = np.zeros(len(fnlist))
    i = 0
    for fn in fnlist:
        goturms = calc_urms(get_mesh_from_file(fn,256,1))
        urms_list[i] = goturms
        i = i+1


    return urms_list
