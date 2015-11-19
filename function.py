import numpy as np
"""
These functions form the basis to solve the 4 FWM equations in Agrawal,
what needs to change is the overlap integrals. As for now they only estimate
1/Aeff

"""

def dbm_to_w(dBm):
    """This function converts a power given in dBm to a power given in W."""
    return 1e-3*10**((dBm)/10.)


def over1(i, k):
    if mode == 'lpo1':
        return (161e-12)**(-1)
    else:
        return (170e-12)**(-1)


def over2(i, j, k, l):
    if mode == 'lpo1':
        return (161e-12)**(-1)
    else:
        return (170e-12)**(-1)


def sums(i, A):
    sum_full = 0
    for j in range(4):
        if j==i:
            break
        sum_full += over1(i,j)*np.abs(A[j])**2
    return sum_full


def integrand(i, A, zz):
        temp = (over1(i, i)*np.abs(A[i])**2+2*sums(i, A))*A[i]
        i_vec = np.array([[0, 1, 2, 3], [1, 0, 2, 3], [2, 3, 0, 1], [3, 2, 0, 1]])
        if 2//(i+1) == 0:
            temp += 2*over2(i_vec[i,0],i_vec[i,0],i_vec[i,0],i_vec[i,0])*\
                np.conjugate(A[i_vec[i,1]])*A[i_vec[i,2]]*A[i_vec[i,3]]*np.exp(-1j*Dk*zz)
        else:
            temp += 2*over2(i_vec[i,0],i_vec[i,0],i_vec[i,0],i_vec[i,0])*\
                np.conjugate(A[i_vec[i,1]])*A[i_vec[i,2]]*A[i_vec[i,3]]*np.exp(1j*Dk*zz)
        temp +=0.5*a*A[i]
        return 1j*n2*omega[i]*temp/c

def system(zz,AB):
    dABdz=np.zeros(4,dtype='complex')
    A0 = AB[0] + AB[1]*1j
    A1 = AB[2] + AB[3]*1j
    A2 = AB[4] + AB[5]*1j
    A3 = AB[6] + AB[7]*1j
    A = [A0, A1, A2, A3]
    dABdz[0] = integrand(0,A,zz)
    dABdz[1] = dABdz[0]#integrand(0,A,zz)
    dABdz[2] = integrand(2,A,zz)
    dABdz[3] = integrand(3,A,zz)
    return dABdz.view(np.float64)

def system2(AB,zz):
    dABdz=np.zeros(4,dtype='complex')
    A0 = AB[0] + AB[1]*1j
    A1 = AB[2] + AB[3]*1j
    A2 = AB[4] + AB[5]*1j
    A3 = AB[6] + AB[7]*1j
    A = [A0, A1, A2, A3]
    dABdz[0] = integrand(0,A,zz)
    dABdz[1] = dABdz[0]#integrand(1,A,zz)
    dABdz[2] = integrand(2,A,zz)
    dABdz[3] = integrand(3,A,zz)
    return dABdz.view(np.float64)