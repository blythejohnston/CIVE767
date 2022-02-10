#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from numpy import sqrt
from matplotlib import pyplot as plt


# In[2]:


#Variables
m = 20.0 #mass
k = 200.0 #stiffness
ksi = 0.0 #damping ratio
wn = sqrt(k/m)
c = ksi*2*sqrt(m*k) # determine damping coef. from damping ratio 
T1=100 #duration of motion evaluation


# In[3]:


#Initial Conditions
F0=0 #no forcing function in free vibration
x0=2 #initial displacement
v0=1 # initial velocity
#x=np.array([x0,v0]) #Velocity,displacement


# In[4]:


#Setup Sampling
deltaT=0.02 #time step
N=T1/deltaT #total number of time points
t=np.arange(0.0,T1,deltaT) #generate time vector
x = []
v = []
X = []
for ti in t:
  x = (x0)*np.cos(wn*ti)+(v0/wn)*np.sin(wn*ti)
  X.append(x)
print(wn)


# In[5]:


plt.plot(t,X)
plt.show()


# In[ ]:




