#!/usr/bin/env python
# coding: utf-8

# ## SDOF Damped under Free Vibration

# ### Import Necessary Packages and Functions

# In[1]:


import numpy as np
from numpy import sqrt
from matplotlib import pyplot as plt


# ### Set Variable Values for System

# In[2]:


m = 50 #mass
k = 51.5 #stiffness
xi = 0.01 #damping ratio - use decimal not percentage
wn = sqrt(k/m)
c = xi*2*sqrt(m*k) # determine damping coef. from damping ratio
wd = wn*sqrt(1.0-xi**2.0)


# ### Establish Initial Conditions

# In[3]:


F0=0 #no forcing function in free vibration
x0=2 #initial displacement
v0=1 # initial velocity


# ### Setup Sampling Protocol

# In[4]:


T1=100 #duration of motion evaluation
deltaT=0.02 #time step
N=T1/deltaT #total number of time points
t=np.arange(0.0,T1,deltaT) #generate time vector
x = []
v = []
X = []
for ti in t:
  x = (np.exp(-xi*wn*ti))*((x0)*np.cos(wd*ti)+((v0+xi*wn*x0)/wd)*np.sin(wd*ti))
  X.append(x)


# ### Plotting the Resulting Displacement Response

# In[5]:


plt.plot(t,X)
plt.show()

