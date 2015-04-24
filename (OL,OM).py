from __future__ import division
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


#Import Data

path = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/SCPUnion2_mu_vs_z.txt'
data = np.genfromtxt(path)


# In[884]:

#Define important arrays
z = data[:,1] # redshift
m = data[:,2] # magnitude
e = data[:,3] # error



## OMEGA M OMEGA DE

# In[ ]:

def mth2(z,omega_m,omega_l,M):  #fixed omega_m
    def dl(z,omega_m,omega_l):
        omega_k = omega_m + omega_l - 1
        f = lambda x: 1/np.sqrt(omega_m*(1+x)**3 + omega_l - omega_k*(1+x)**2)
        d = integrate.romberg(f, 0 , z, rtol = 0.001, show = False)*(1+z)
        return d 
    return 5*log10(dl(z,omega_m,omega_l)) + M 


# In[827]:

#fit script M again

#Create Arrays for low z values

thres = .07
total = zip(z,m,e)
new = filter(lambda total: total[0]<thres,total)
znew, mnew, enew = zip(*new)

#Create array for M
Mi = 43
Mf = 43.5
dM = .0001
Mmax = np.arange(Mi,Mf,dM)

chi2 = np.zeros(len(Mmax))

j = 0
for a in Mmax: # runs through differen M values
    for i in range(len(znew)):        # runs through different supernova
        chi2[j] += (mnew[i]-mth2(znew[i],.3,.7,a))**2/enew[i]**2
    #print  chi2[j]
    j +=1

L = np.exp(-chi2/2)
normL5 = L/np.max(L)

#find max likelihood for given M

i = np.unravel_index(normL5.argmax(),shape(normL5))
M = Mi + i[0]*dM

#calculate standard deviation
for i in range(len(Mmax)):
    if normL3[i] > .35:
        std5 = abs(Mi+i*dM-M)
        break

print "Best fit M = ",M,"+/- ",std5


# In[720]:

M = 43.165 # found from above fitting
domega = .002

omega_l_i = 0
omega_l_f = 1
omega_l_max = np.arange(omega_l_i,omega_l_f,domega)

omega_m_i = 0
omega_m_f = 1
omega_m_max = np.arange(omega_m_i,omega_m_f,domega)

chi2 = np.zeros([len(omega_m_max),len(omega_l_max)])

i = 0
for omega_m in omega_m_max: 
    print omega_m
    j = 0
    for omega_l in omega_l_max:
        for k in range(len(z)):        # runs through different supernova
            chi2[i,j] += ((m[k]-mth2(z[k],omega_m,omega_l,M))/e[k])**2
        j += 1
    i += 1

L = np.exp(-chi2/2)
normL6 = L/np.max(L)

i,j = np.unravel_index(normL6.argmax(),normL6.shape)
OL = omega_l_i + j*domega
OM = omega_m_i + i*domega

print "Best fit omega_m = ",OM
print "Best fit omega_l = ",OL

test2 = np.zeros(shape(normL6))
for i in range(shape(normL6)[0]):
    for j in range(shape(normL6)[1]):
        if normL6[i,j] > .32:
            test2[i,j] = 1
        elif normL6[i,j] > .05 and normL6[i,j] <.32:
            test2[i,j] = 1/2
        else:
            test2[i,j] = 0

numpy.savetxt("OMmax2.csv",omega_m_max,delimiter = ",")
numpy.savetxt("OLmax.csv",omega_l_max,delimiter = ",")
numpy.savetxt("OLvOM.csv", normL5, delimiter=",")
numpy.savetxt("OLvOMtest.csv", test2, delimiter=",")


# In[891]:

###heres another plot comparing the theoretical magnitude in different cosmologies. 

zmax = np.arange(.01,1.5,.01)
val = np.zeros([3,len(zmax)])

OM = 0.264
OL = .708
w = -1
M = 43.164

for i in range(len(zmax)):
    val[0,i]=mth2(zmax[i],OM,OL,M)
    val[1,i]=mth2(zmax[i],1,0,M)
    val[2,i]=mth2(zmax[i],0,1,M)
    
numpy.savetxt("magvred2.csv",val,delimiter = ",")

