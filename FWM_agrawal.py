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
def mw_to_dbm(W):
    """This function converts a power given in mW to a power given in dBm."""
    W *=1e-3    
    return 10.*np.log10(W)

def sums(i, A):
    sum_full = 0
    for j in range(4):
        if i != j:           
           sum_full += np.abs(np.sign(i-j))*over1()*np.abs(A[j])**2
    return sum_full


def integrand(i, A, zz):
        temp = (over1()*np.abs(A[i])**2+2*sums(i, A))*A[i]
        i_vec = np.array([[1, 2, 3], [0, 2, 3], [3, 0, 1], [2, 0, 1]])
                
        if 2//(i+1) == 0:
            #print 'signal or idler', i            
            temp += 2*over1()*\
                np.conjugate(A[i_vec[i,0]])*A[i_vec[i,1]]*A[i_vec[i,2]]*np.exp(-1j*Dk*zz)
        else:
            #print 'pump', i      
            temp += 2*over1()*\
                np.conjugate(A[i_vec[i,0]])*A[i_vec[i,1]]*A[i_vec[i,2]]*np.exp(1j*Dk*zz)
        #temp -=0.5*a*A[i]
        return 1j*n2*omega[i]*temp/c

    
def system(zz,AB):
    try:
        temp = AB[0]
    except:
        temp = AB
        AB = np.copy(zz)        
        zz = np.copy(temp)
        del temp
        pass
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


a = 0 # attenuation of the fibre
n2 = 2.35e-26 # nonlinear coefficient
num_steps = 10000 # number of steps for the integrator
d = 1e3 # propagation distance
mode='lpo1' # mode that is simmulated
lamda_c = 1.5508e-6
zmin =0 #starting point of the fibre
zmax = d #end point of interest

if mode == 'lpo1':
    def over1():
        return (161e-12)**(-1)
else:
    def over1():
        return (170e-12)**(-1)

if mode == 'lpo1':
    D = 19.9*1e-6
    beta2 = -lamda_c**2*D/(2*pi*c) 
else:
    D = 22*1e-6
    beta2 = -lamda_c**2*D/(2*pi*c)

P_vec = np.arange(23,32,2)
lamp_vec = np.array([1.549])

lams_vec = np.linspace(1.548,1.552,1000)
lami_vec = np.zeros(len(lams_vec))

AB_final = np.zeros([4,len(P_vec),len(lams_vec)],dtype='complex')
lamp = lamp_vec[0]
lami = lami_vec[0]

z = np.linspace(zmin,zmax,num_steps) #linearly spaced propagation vector

dz = z[1] - z[0]
for i,Po in enumerate(P_vec):
    P = dbm_to_w(Po)  
    AB0 = np.array([P,P, P/100, P/1000], dtype='complex')
            
    for j,lams in enumerate(lams_vec):            
        lam_vec = np.array([lamp,lamp,lams,lami])
        lam_vec *=1e-6
        omega= 2*pi*c/lam_vec
        omega[3] = 2*omega[0]-omega[2]
        lami_vec[j] = 2*pi*c/omega[3]
        f = omega/(2*pi)            
        
                     
        Dk = np.abs(0.5*beta2*(omega[2]**2+omega[3]**2-2*omega[0]**2))
        r = ode(system).set_integrator('dop853')
        r.set_initial_value(AB0.view(np.float64),np.float64(zmin))
        AB_vec = np.zeros([num_steps,4],dtype='complex')
        count = 1
        AB_vec[0,:] = AB0
        
        while count < num_steps and r.successful():
            r.integrate(r.t+dz)
            count+=1
        if r.successful():
            AB_final[:,i,j] = r.y.view(np.complex128)
        else:
            print("object integration failed, trying the lsoda")          
            AB1 = np.zeros([len(z), 3])
            AB1 = odeint(system, AB0.view(np.float64), z,atol=None,rtol=None,full_output = 1)
            if AB1[1]['message']  == 'Integration successful.':
                AB_final[:,i,j] =  np.shape((AB1[0][-1::].view(np.complex)))
            else:
                sys.exit('both integrations failed')
print np.abs(AB0)**2 - np.abs(AB_final)**2
sys.exit()
for i in P_vec:
    plt.plot(lams_vec,mw_to_dbm(np.abs(AB_final[2,0,:])**2),label=i )
plt.plot(lams_vec,mw_to_dbm(np.abs(AB0[2])**2),'o',label='before')
plt.xlabel(r'$\lambda (\mu m)$')
plt.title('Signal')
plt.ylabel(r'$P(W)$')
plt.legend()
plt.show()
#lami_vec = 2*pi*c/omega[3]
plt.plot(lami_vec,mw_to_dbm(np.abs(AB_final[3,0,:])**2),'o',label='after')
#plt.plot(lami_vec,mw_to_dbm(np.abs(AB0[3])**2),'o',label='before')
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