
# coding: utf-8

# In[902]:

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


# In[910]:

#Import Data

path = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/SCPUnion2_mu_vs_z.txt'
data = np.genfromtxt(path)


# In[904]:

#Define important arrays
z = data[:,1] # redshift
m = data[:,2] # magnitude
e = data[:,3] # error


## OMEGA M AND W

# In[905]:

#define a function that calculates the theoretical magnitude based on cosmological parameters

def mth(z,omega_m,M,w):  #fixed omega_m
    def dl(z,omega_m,w):
        f = lambda x: 1/np.sqrt(omega_m*(1+x)**3 + (1-omega_m)*(1+x)**(3*(1+w)))
        d = integrate.romberg(f, 0 , z, rtol = 0.001, show = False)*(1+z)
        return d
    return 5*log10(dl(z,omega_m,w)) + M


#### Finding Script M

# In[906]:

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

print "Best fit M = ",M,"+/- ",std1


#Save data to csv for plotting later
LvW = zip(Mmax,normL1)
numpy.savetxt("LvM.csv", LvW, delimiter=",")


####  Fixed w = -1, varying Omega_m  

# In[907]:

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

print "Best fit Omega_M = ",OM, "+/- ",std2

LvW = zip(omega_m_max,normL2)
numpy.savetxt("LvOM.csv", LvW, delimiter=",")


####  Fixed Omega_m = .33, varying w  

# In[908]:

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

print "Expected value of w = ",W, "+/- ", std3

LvW = zip(wmax,normL3)
numpy.savetxt("LvW.csv", LvW, delimiter=",")


##### # Here I realized that the value of script M that I had fit before is close to the best fit value, but not quite, I had fit Script M while simultaneously fitting the other parameters (omega_m and w). The gridding for M is quite fine while the gridding for the other two parameters is large. 

# In[ ]:

#Finding script M again, with all three parameters


# Create arrays
dw = .01
wi = -1.4
wf = -.6 + dw
wmax = np.arange(wi,wf,dw)

domega = .005
omega_m_i = .25
omega_m_f = .4 + domega
omega_m_max = np.arange(omega_m_i,omega_m_f,domega)

Mi = 43.1
Mf = 43.2
dM = .005
Mmax = np.arange(Mi,Mf,dM)

chi2 = np.zeros([len(omega_m_max),len(wmax),len(Mmax)])


#run the chi2 loop
n = 0
for a in Mmax:
    print "M = ",a
    i = 0
    for omega in omega_m_max: 
        print "Omega_m = ",omega
        j = 0
        for w in wmax:
            #print "w = ",w
            for k in range(len(z)):        # runs through different supernova
                chi2[i,j,n] += ((m[k]-mth(z[k],omega,a,w))/e[k])**2
            j += 1
        i += 1
    n += 1

L = np.exp(-chi2/2)
normL4 = L/np.max(L)

i,j,n = np.unravel_index(normL4.argmax(),normL4.shape)
W = wi + j*dw
OM = omega_m_i + i*domega
M = Mi + n*dM

print "Best fit omega_m = ",OM
print "Best fit w = ",W
print "Best fit M = ",M


##### Thus, my new, better fitting Script M value is 43.15 rather than 43.2. This makes quite a difference. I accept this value as best fit and then vary Omega_m and w in one plane, rather than iterating over every M value. This was because gridding over all three parameters would simply take too long. 

#### Varying w and omega_m with fixed M

# In[ ]:

### Takes FOREVER to run (6 hours) ###

M = 43.155 # found by fitting all three parameters

dw = .001
wi = -1.4
wf = -.6 + dw
wmax = np.arange(wi,wf,dw)

domega = .001
omega_m_i = 0
omega_m_f = .5 + domega

omega_m_max = np.arange(omega_m_i,omega_m_f,domega)

chi2 = np.zeros([len(omega_m_max),len(wmax)])

i = 0
for omega in omega_m_max: 
    print omega
    j = 0
    for w in wmax:
        for k in range(len(z)):        # runs through different supernova
            chi2[i,j] += ((m[k]-mth(z[k],omega,M,w))/e[k])**2
        j += 1
    i += 1

L = np.exp(-chi2/2)
normL4 = L/np.max(L)

i,j = np.unravel_index(normL4.argmax(),normL4.shape)
W = wi + j*dw
OM = omega_m_i + i*domega

print "Best fit omega_m = ",OM
print "Expected value of w = ",W


test1 = np.zeros(shape(normL4))
for i in range(shape(normL4)[0]):
    for j in range(shape(normL4)[1]):
        if normL4[i,j] > .32:
            test1[i,j] = 1
        elif normL4[i,j] > .05 and normL4[i,j] <.32:
            test1[i,j] = 1/2
        else:
            test1[i,j] = 0

numpy.savetxt("OMmax1.csv",omega_m_max,delimiter = ",")
numpy.savetxt("Wmax1.csv",wmax,delimiter = ",")
numpy.savetxt("WvOM.csv", normL4, delimiter=",")
numpy.savetxt("WvOMtest.csv", test1, delimiter=",")


##### Here's a plot with three different magnitude curves. 

# In[ ]:

zmax = np.arange(.01,1.5,.01)
val = np.zeros([3,len(zmax)])

w = -1.04
for i in range(len(zmax)):
    val[0,i]=mth(zmax[i],OM,M,w)
    val[1,i]=mth(zmax[i],0,M,w)
    val[2,i]=mth(zmax[i],1,M,w)

numpy.savetxt("magvred.csv",val,delimiter = ",")
numpy.savetxt("zmax.csv",zmax,delimiter = ",")


## OMEGA M OMEGA DE

# In[909]:

def mth2(z,omega_m,omega_l,M):  #fixed omega_m
    def dl(z,omega_m,omega_l):
        omega_k = omega_m + omega_l - 1
        f = lambda x: 1/np.sqrt(omega_m*(1+x)**3 + omega_l - omega_k*(1+x)**2)
        d = integrate.romberg(f, 0 , z, rtol = 0.001, show = False)*(1+z)
        return d 
    return 5*log10(dl(z,omega_m,omega_l)) + M 


# In[ ]:

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
##
print "Best fit M = ",M,"+/- ",std5


# In[ ]:

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


# In[ ]:

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


# In[ ]:



