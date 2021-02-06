#!/usr/bin/env python
# coding: utf-8

# # SMIB system as in Milano's book example 8.1

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from IPython.core.display import HTML
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")
plt.ion()


# In[2]:


get_ipython().run_line_magic('matplotlib', 'widget')


# In[3]:


from pydae import ssa
from smib_milano_ex8p1_4ord import smib_milano_ex8p1_4ord_class


# ### Instantiate system

# In[4]:


smib = smib_milano_ex8p1_4ord_class()
smib.t_end = 15.0
smib.Dt = 0.01
smib.decimation =1
smib.update()


# ### Initialize the system (backward and foreward)

# In[5]:


events=[{'p_t':0.1,'X_line':0.01}]
smib.initialize(events,xy0='xy_0.json')
smib.save_0()

smib.report_x()  # obtained dynamic states
smib.report_y()  # obtained algebraic states
smib.report_z()  # obtained outputs
smib.report_u()  # obtained algebraic states (theta is both state and output; f_x is both input and output)
smib.report_params()  # considered parameters


# ### Small signal analisys

# In[6]:


ssa.eval_A(smib)              # method to linealized the system and to compute matrix A
eig_df=ssa.damp_report(smib)  # method to create a pandas.DataFrame after computing eigenvalues for A
eig_df


# In[7]:


ssa.participation(smib).abs().round(2)  


# ### Simulation
# 
# A time simulation can be performed using the method `simulate`:

# In[8]:


smib.initialize([{'p_t':0.7, 'q_t':0.2}],xy0='xy_0.json')
smib.simulate([{'t_end':1},   # compute initial condition with defined P and Q and run until t=1s
               {'t_end':20, 'v_f':2.5}],           # compute initial condition with defined P and Q and run until t=1s   
                'prev');                              


# In[9]:


# plotting the results with matplolib:
plt.close('all')
fig, axes = plt.subplots(nrows=1,ncols=1, figsize=(6, 3), dpi=100)

axes.plot(smib.T, smib.get_values('p_t') , label=f'$p_t$')
axes.plot(smib.T, smib.get_values('q_t') , label=f'$q_t$')
axes.grid()
axes.set_ylabel('Powers (p.u.)')
axes.set_xlabel('Time (s)')
axes.legend()
fig.tight_layout()


# In[ ]:





# In[10]:


def simulate_fault(system,duration,X_fault=1e-4,N_steps=500):
    Dt = smib.struct[0].Dt
    Dt_recovery = 10e-6
    t_0 = smib.struct[0].t
    system.run([{'t_end':t_0+duration,'B_shunt':1/X_fault}])
    it = 0.0
    for x in np.linspace(0,1.0,N_steps):
        #B_shunt = 1/X_fault - x**2/X_fault
        B_shunt = 1/X_fault - x/X_fault
        
        B_shunt = 1/X_fault - x**0.5/X_fault
        system.run([{'Dt':Dt_recovery/2,'t_end':t_0+duration+it*Dt_recovery,'B_shunt':B_shunt}])
        it+=1.0
    system.run([{'Dt':Dt_recovery/2,'t_end':t_0+duration+(it+1)*Dt_recovery,'B_shunt':0.0}])
    system.run([{'Dt':Dt,'t_end':t_0+duration+Dt,'B_shunt':0.0}])
    
    
smib = smib_milano_ex8p1_4ord_class()

smib.initialize([{'p_t':0.8, 'q_t':0.4, 'D':0.0}],1)
smib.run([{'t_end':1}])  
simulate_fault(smib,0.05)
smib.run([{'t_end':20}])

smib.post();


# In[11]:


# plotting the results with matplolib:
plt.close('all')
fig, axes = plt.subplots(nrows=3,ncols=1, figsize=(6, 6), dpi=100, sharex=True)

axes[0].plot(smib.T, smib.get_values('v_1') , label=f'$v_1$')

axes[1].plot(smib.T, np.rad2deg(smib.get_values('delta')) , label=f'$\delta$')

axes[2].plot(smib.T, smib.get_values('omega') , label=f'$\omega$')

axes[0].set_ylabel('Powers (p.u.)')
axes[1].set_xlabel('Time (s)')


for ax in axes:
    ax.legend()
    ax.grid()
    
fig.tight_layout()


# In[12]:


# plotting the results with matplolib:
plt.close('all')
fig, axes = plt.subplots(nrows=1,ncols=1, figsize=(6, 3), dpi=100)

axes.plot(np.rad2deg(smib.get_values('delta')), smib.get_values('p_m') , label=f'$p_m$')
axes.plot(np.rad2deg(smib.get_values('delta')), smib.get_values('p_t') , label=f'$p_t$')

axes.grid()
axes.set_ylabel('Powers (p.u.)')
axes.set_xlabel('Time (s)')
axes.legend()
fig.tight_layout()


# In[13]:


smib.decimation


# In[14]:


smib.struct[0]['iters'][smib.struct[0]['it']]


# In[139]:





# In[157]:


X_fault = 1e-3
x = np.linspace(0,1.0,200) 
B_shunt = 1/X_fault - x**0.5/X_fault

fig, axes = plt.subplots(nrows=1,ncols=1, figsize=(6, 3), dpi=100)

axes.plot(x, B_shunt)


# In[ ]:





# In[74]:


import xml.etree.ElementTree


# In[ ]:




