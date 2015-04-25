
# coding: utf-8

# In[196]:

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
import matplotlib 


# In[197]:

#Import Data

#from numpy import genfromtxt
#my_data = genfromtxt('my_file.csv', delimiter=',')

path1 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/SCPUnion2_mu_vs_z.txt'

path2 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/LvM.csv'
path3 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/LvOM.csv'
path4 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/LvW.csv'

path5 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/OMmax1.csv'
path6 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/Wmax1.csv'
path7 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/WvOM.csv'
path8 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/WvOMtest.csv'

path9 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/OLmax.csv'
path10 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/OMmax2.csv'
path11 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/OLvOM.csv'
path12 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/OLvOMtest.csv'

path13 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/magvred.csv'
path14 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/magvred2.csv'
path15 = '/Users/Manoj/Dropbox/Umich/W2015/P457/Project/Physics-457-Supernova-Project/My Data/zmax.csv'

data = np.genfromtxt(path1)

LvM = np.genfromtxt(path2,delimiter = ',')
LvOM = np.genfromtxt(path3,delimiter = ',')
LvW = np.genfromtxt(path4,delimiter = ',')

OMmax1 = np.genfromtxt(path5,delimiter = ',')
Wmax1 = np.genfromtxt(path6,delimiter = ',')
WvOM = np.genfromtxt(path7,delimiter = ',')
WvOMtest = np.genfromtxt(path8,delimiter = ',')

OLmax = np.genfromtxt(path9,delimiter = ',')
OMmax2 = np.genfromtxt(path10,delimiter = ',')
OLvOM = np.genfromtxt(path11,delimiter = ',')
OLvOMtest = np.genfromtxt(path12,delimiter = ',')

magvred = np.genfromtxt(path13,delimiter = ',')
magvred2 = np.genfromtxt(path14,delimiter = ',')
zmax = np.genfromtxt(path15,delimiter = ',')



# In[198]:

#Define important arrays

z = data[:,1] # redshift
m = data[:,2] # magnitude
e = data[:,3] # error


### Figure 1

# In[199]:

figure(num=None, figsize=(10, 6), dpi=280, facecolor='w', edgecolor='k')

X,Y = zip(*LvM)
plt.plot(X,Y,".")
plt.title("Likelihood ($\mathcal{L}$) vs. $w$",fontsize = 35)
plt.xlabel("$\mathcal{M}$",fontsize = 30)
plt.ylabel("Likelihood ($\mathcal{L}$)",fontsize = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()


### Figure 2a

# In[200]:

figure(num=None, figsize=(10, 6), dpi=280, facecolor='w', edgecolor='k')

X,Y = zip(*LvOM)
plt.plot(X,Y,".")
plt.title("Likelihood ($\mathcal{L}$) vs. $\Omega_M$",fontsize = 35)
plt.xlabel("$\Omega_M$",fontsize = 30)
plt.ylabel("Likelihood ($\mathcal{L}$)",fontsize = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()


### Figure 2b

# In[201]:

figure(num=None, figsize=(10, 6), dpi=280, facecolor='w', edgecolor='k')

X,Y = zip(*LvW)
plt.plot(X,Y,".")
plt.title("Likelihood ($\mathcal{L}$) vs. $w$",fontsize = 35)
plt.xlabel("$w$",fontsize = 30)
plt.ylabel("Likelihood ($\mathcal{L}$)",fontsize = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()


### Figure 3

# In[202]:

X, Y = np.meshgrid(Wmax1,OMmax1)


figure(num=None, figsize=(10, 8), dpi=580, facecolor='w', edgecolor='k')

plt.contourf(Y,X,WvOMtest[:500,:], 32, rstride=1, cstride=1, cmap=cm.Blues)
#plt.text(0.5, 0.5,'$w$ = -1.04, $\Omega_m$ = 0.284',fontsize = 20, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.title("Likelihood Surface ($\mathcal{L}$) for $w$ vs. $\Omega_m$",fontsize = 30)
plt.ylabel("$w$",fontsize = 30)
plt.xlabel("$\Omega_M$",fontsize = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()


### Figure 4

# In[203]:

figure(num=None, figsize=(10, 6), dpi=280, facecolor='w', edgecolor='k')
plt.errorbar(z,m,yerr=e,fmt='.',color = 'b')


plt.plot(zmax,magvred[0,:],label = "Omega_m = .33")
plt.plot(zmax,magvred[1,:],label = "Omega_m = 0")
plt.plot(zmax,magvred[2,:],label = "Omega_m = 1")
plt.title("Magnitude vs. Redshift ($w,\Omega_m$)", fontsize = 30)
plt.xlabel("Redshift",fontsize = 20)
plt.ylabel("Magnitude",fontsize = 20)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.legend(loc=4,prop={'size':15},fancybox = True)
plt.show()


### Figure 5

# In[204]:

X, Y = np.meshgrid(OLmax,OMmax2)


figure(num=None, figsize=(10, 8), dpi=580, facecolor='w', edgecolor='k')

plt.contourf(Y,X,OLvOMtest, 32, rstride=1, cstride=1, cmap=cm.Blues)
#plt.text(0.5, 0.5,'$w$ = -1.04, $\Omega_m$ = 0.284',fontsize = 20, horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)
plt.title("Likelihood Surface ($\mathcal{L}$) for $\Omega_\Lambda$ vs. $\Omega_m$",fontsize = 30)
plt.ylabel("$\Omega_\Lambda$",fontsize = 25)
plt.xlabel("$\Omega_M$",fontsize = 25)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()


### Figure 6

# In[205]:

figure(num=None, figsize=(10, 6), dpi=280, facecolor='w', edgecolor='k')
plt.errorbar(z,m,yerr=e,fmt='.',color = 'b')


plt.plot(zmax,magvred2[0,:],label = "Best fit")
plt.plot(zmax,magvred2[1,:],label = "$\Omega_m$ = 1")
plt.plot(zmax,magvred2[2,:],label = "$\Omega_{DE}$ = 1")
plt.title("Magnitude vs. Redshift ($\Omega_{DE},\Omega_m$)", fontsize = 30)
plt.xlabel("Redshift, $z$",fontsize = 20)
plt.ylabel("Magnitude",fontsize = 20)
plt.legend(loc=4,prop={'size':15},fancybox = True)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.show()

