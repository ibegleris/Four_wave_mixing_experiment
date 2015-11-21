from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from scipy.constants import pi,c
from scipy.integrate import odeint, ode
from scipy.io import loadmat
import sys


def dbm_to_w(dBm):
    """This function converts a power given in dBm to a power given in W."""
    return 1e-3*10**((dBm)/10.)


def over1():
    """Finds the overlap integrals, since only one mode is simmulated at the time
    the ovelap is equal to the inverse of A_eff"""
    if mode == 'lpo1':
        return (161e-12)**(-1)
    else:
        return (170e-12)**(-1)

def sums(i, A):
    sum_full = 0.
    for j in range(3):
        if j != i:
            sum_full +=np.abs(A[j])**2
    return sum_full


def integrand(i, A, zz):
        temp = (np.abs(A[i])**2+2*sums(i, A))*A[i]
            
        if 2//(i+1) == 2:
            temp += 2*np.conjugate(A[i])*A[i+1]*\
                    A[i+2]*np.exp(1j*Dk*zz)
        else:
            temp +=2*np.conjugate(A[i])*A[0]**2\
                    *np.exp(-1j*Dk*zz)
        temp -=0.5*a*A[i]
        gama = over1()*1j*n2*omega[0]/c    
        return temp*gama

def system(zz,AB):
    #print type(zz)
    try:
        temp = AB[0]
    except:
        temp = AB
        AB = np.copy(zz)        
        zz = np.copy(temp)
        del temp
        pass
    dABdz=np.zeros(3,dtype='complex')
    A0 = AB[0] + AB[1]*1j
    A1 = AB[2] + AB[3]*1j
    A2 = AB[4] + AB[5]*1j
    A = [A0, A1, A2]
    dABdz[0] = integrand(0,A,zz)
    dABdz[1] = integrand(1,A,zz)
    dABdz[2] = integrand(2,A,zz)
    return dABdz.view(np.float64)


#P = 18 #propagating power in dbm
#P = dbm_to_w(P) # convert to W
a = 0 # attenuation of the fibre
n2 = 2e-20 # nonlinear coefficient
num_steps = 1000000 # number of steps for the integrator
d = 1e3 # propagation distance
mode='lpo1' # mode that is simmulated
if mode == 'lpo1':
    D = 19.9*1e-6
else:
    D = 14.4*1e-6
P_vec = np.arange(23,33,2)
lamp_vec = np.array([1.5508])
lams_vec = np.array([1.5509])
lami_vec = np.array([1.560])
AB_final = np.zeros([3,len(lamp_vec),len(lams_vec),len(lami_vec)],dtype='complex')
for Po in P_vec:
    P = dbm_to_w(Po)  
    for i,lamp in enumerate(lamp_vec):
        for j,lams in enumerate(lams_vec):    
            for k,lami in enumerate(lami_vec):            
                
                mat = loadmat('group_velocities_inv.mat')
                if mode == 'lpo1':
                    dk = mat['invg1']
                else:
                    dk = mat['invg2']
                AB0 = np.array([P, P/10, 0], dtype='complex')
                Dk=1 # Phase Matching
                lam_vec = np.array([lamp,lams,lami])
                lam_vec *=1e-6
                omega= 2*pi*c/lam_vec
                omega[2] = 2*omega[0]-omega[1]
                f = omega/(2*pi)            
                Dk= -8*pi*D*f[0]**2*(f[1]-f[0])/c
                sys.exit()
                zmin =0 #starting point of the fibre
                zmax = d #end point of interest
                z = np.linspace(zmin,zmax,num_steps) #linearly spaced propagation vector
                dz = z[1] - z[0]
                
                r = ode(system).set_integrator('dop853')
                r.set_initial_value(AB0.view(np.float64),np.float64(zmin))
                AB_vec = np.zeros([num_steps,3],dtype='complex')
                count = 1
                AB_vec[0,:] = AB0
                while count < num_steps and r.successful():
                    r.integrate(r.t+dz)
                    #AB_vec[count,:] = r.y.view(np.complex128)
                    count+=1
                
                if r.successful() == False:
                    print("object integration failed, trying the lsoda")          
                    AB1 = np.zeros([len(z), 3])
                    AB1 = odeint(system, AB0.view(np.float64), z,atol=None,rtol=None,full_output = 1)
                    if AB1[1]['imxer'] >0:
                        AB_final[:,i,j,k] = AB1[0][-1::,:].view(np.complex128)
                    else:
                        sys.exit('both integrations failed')
                else:
                    AB_final[:,i,j,k] = r.y.view(np.complex128)


plt.plot(lamp_vec,np.abs(AB_final[0,:,0,0])**2,'o',label='after')
plt.plot(lamp_vec,np.abs(AB0[0])**2,'o',label='before')
plt.xlabel(r'$\lambda (\mu m)$')
plt.title('Pump')
plt.ylabel(r'$P(W)$')
plt.legend()
plt.show()

plt.plot(lams_vec,np.abs(AB_final[1,:,0,0])**2,'o',label='after')
plt.plot(lams_vec,np.abs(AB0[1])**2,'o',label='before')
plt.xlabel(r'$\lambda (\mu m)$')
plt.title('Signal')
plt.ylabel(r'$P(W)$')
plt.legend()
plt.show()

plt.plot(lami_vec,np.abs(AB_final[2,:,0,0])**2,'o',label='after')
plt.plot(lami_vec,np.abs(AB0[2])**2,'o',label='before')
plt.xlabel(r'$\lambda (\mu m)$')
plt.title('idler')
plt.ylabel(r'$P(W)$')
plt.legend()
plt.show()



"""
plt.plot(z,np.abs(AB1[:,0])**2,label='pump')
plt.plot(z,np.abs(AB1[:,1])**2,label='probe')
plt.plot(z,np.abs(AB1[:,2])**2,label='idler')
plt.legend()

plt.title('First int')
plt.show()


plt.plot(z,np.abs(AB_vec[:,0])**2,label='pump')
plt.plot(z,np.abs(AB_vec[:,1])**2,label='probe')
plt.plot(z,np.abs(AB_vec[:,2])**2,label='idler')
plt.legend()
plt.title('Second int')
plt.show()
"""