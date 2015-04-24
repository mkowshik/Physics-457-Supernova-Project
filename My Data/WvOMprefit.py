
from __future__ import division
from __future__ import print_function
from scipy import integrate
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import math
from math import log10,exp
import sys
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm


# In[883]:

#Import Data

path = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/SCPUnion2_mu_vs_z.txt'
data = np.genfromtxt(path)


# In[884]:

#Define important arrays
z = data[:,1] # redshift
m = data[:,2] # magnitude
e = data[:,3] # error


## OMEGA M AND W

# In[892]:

#define a function that calculates the theoretical magnitude based on cosmological parameters

def mth(z,omega_m,M,w):  #fixed omega_m
    def dl(z,omega_m,w):
        f = lambda x: 1/np.sqrt(omega_m*(1+x)**3 + (1-omega_m)*(1+x)**(3*(1+w)))
        d = integrate.romberg(f, 0 , z, rtol = 0.001, show = False)*(1+z)
        return d
    return 5*log10(dl(z,omega_m,w)) + M


#### Finding Script M

# In[890]:

#Create Arrays for low z values where all models converge

thres = .07  #low z threshold

#splice the data down
total = zip(z,m,e)
new = filter(lambda total: total[0]<thres,total)
znew, mnew, enew = zip(*new)

#Create array for M
Mi = 43
Mf = 43.4
dM = .001
Mmax = np.arange(Mi,Mf,dM)

chi2 = np.zeros(len(Mmax))


#Calculate Chi^2
j = 0
for a in Mmax: # runs through differen M values
    for i in range(len(znew)):        # runs through different supernova
        chi2[j] += (mnew[i]-mth(znew[i],1,a,-1))**2/enew[i]**2
    j +=1

L = np.exp(-chi2/2) #calculate likelihood
normL1 = L/np.max(L) #normalize it

#find max likelihood for given M (find the right index)

i = np.unravel_index(normL1.argmax(),shape(normL1)) 
M = Mi + i[0]*dM

#calculate standard deviation
for i in range(len(Mmax)):
    if normL1[i] > .35:
        std1 = abs(Mi+i*dM-M)
        break

print ("Best fit M = ",M,"+/- ",std1)


#Save data to csv for plotting later
LvW = zip(Mmax,normL1)
numpy.savetxt("LvM.csv", LvW, delimiter=",")


####  Fixed w = -1, varying Omega_m  

# In[889]:

omega_m_i = .2
omega_m_f = .5
domega = .001
omega_m_max = np.arange(omega_m_i,omega_m_f+domega,domega)
    
chi2 = 0*omega_m_max

j = 0
for k in omega_m_max:  # runs through p_j, different cosmo
    for i in range(len(z)):        # runs through different supernova
        chi2[j] += ((m[i]-mth(z[i],k,M,-1))/e[i])**2
    j += 1 

L = np.exp(-chi2/2)
normL2 = L/np.max(L)

#Find max likelihood for given Omega_M

i = np.unravel_index(normL2.argmax(),shape(normL2))
OM = omega_m_i + i[0]*domega

#calculate standard deviation
for i in range(len(omega_m_max)):
    if normL2[i] > .35:
        std2 = abs(omega_m_i+i*domega-OM)
        break

print ("Best fit Omega_M = ",OM, "+/- ",std2)

LvW = zip(omega_m_max,normL2)
numpy.savetxt("LvOM.csv", LvW, delimiter=",")


####  Fixed Omega_m = .33, varying w  

# In[888]:

wi = -3
wf = 3
dw = .01
wmax = np.arange(wi,wf+dw,dw)
    
chi2 = 0*wmax

j = 0
for w in wmax:  # runs through p_j, different cosmo
    for i in range(len(z)):        # runs through different supernova
        chi2[j] += ((m[i]-mth(z[i],OM,M,w))/e[i])**2
    j += 1 

L = np.exp(-chi2/2)
normL3 = L/np.max(L)

#Find most likely w value
i = np.unravel_index(normL3.argmax(),shape(normL3))
W = wi + i[0]*dw

#calculate standard deviation
for i in range(len(wmax)):
    if normL3[i] > .35:
        std3 = abs(wi+i*dw-W)
        break

print ("Expected value of w = ",W, "+/- ", std3)

LvW = zip(wmax,normL3)
numpy.savetxt("LvW.csv", LvW, delimiter=",")
