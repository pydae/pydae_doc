#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'widget')


# In[2]:


import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
import pydae.ssa as ssa
import pydae.grid_tools as gt
from pydae.tools import get_i
import json
import time


# In[3]:


from grid_4bus4wire import grid_4bus4wire_class


# In[4]:


grid = grid_4bus4wire_class()
grid.Dt = 10e-3
grid.decimation = 1
grid.update()


# ## Initialization
# 
# Here the initialization will use the initial guesses from the file `xy_0_dict.json`

# In[5]:


grid.initialize([{}],xy0='xy_0_dict.json')


# The results can be obtained by reporting the algebraic variables:

# In[6]:


grid.report_y()


# As can be observed to get results can be messy by just considering the previous report. 
# Therefore tools can be used from the module `pydae.grid_tools` imported as `gt` here.

# In[7]:


from importlib import reload  
reload(gt)
gt.get_voltage(grid,'B2dc',output='v_an')    # bus B2 voltage 


# In[6]:




events = [{ # CTRL4-3-0
           't_end':0.0, 
                     },
          {'t_end':1.0}, 
          {'t_end':6.0},
          {'t_end':15.0}
          ]

with open('xy_0_dict.json') as fobj:
    data = json.loads(fobj.read().replace("'",'"'))
    
for item in data:
    syst.xy_prev[syst.y_ini_list.index(item)] = data[item]


# In[11]:


V_b = 400/np.sqrt(3)
gt.set_voltages(syst,'B1',[V_b*1,V_b*1,V_b*1.0],0)
gt.set_voltages(syst,'B4',[V_b,V_b,V_b*1],0)
syst.initialization_tol = 1e-12

syst.initialize([{}],xy0='xy_0_dict.json')
t=0


# In[19]:


t+=0.1
gt.set_voltages(syst,'B1',[V_b*1,V_b*1,V_b*1.005],0)
gt.set_voltages(syst,'B4',[V_b,V_b,V_b*0.995],0)
syst.run([{'t_end':t}])

grid_1 = gt.grid(syst)
grid_1.dae2vi()
grid_1.get_v()
grid_1.get_i()

phases = ['a','b','c','n']
print(f'B1:       --- B4:          ')
for it in range(4):
    I_B1 = float(np.abs(grid_1.I_lines[it]))
    I_B4 = float(np.abs(grid_1.I_lines[it+8]))

    print(f'I_{phases[it]} = {I_B1:3.0f} --- I_{phases[it]} = {I_B4:3.0f}')
    

grid_1.bokeh_tools()

gt.plot_results(grid_1)


# In[11]:


syst = grid_4bus4wire_class()
syst.Dt = 10e-3
syst.decimation = 1
syst.update()

Δt = 100.0e-3
Times = np.arange(0,6.0,Δt)

gt.set_voltages(syst,'B1',[V_b*1,V_b*1,V_b*1.0],0)
gt.set_voltages(syst,'B4',[V_b,V_b,V_b*1],0)

syst.initialize([{}],
               'xy_0_dict.json');

N_times = len(Times)

I_B1 = np.zeros((N_times,4))
I_B4 = np.zeros((N_times,4))

xi_V_B1 = np.zeros((N_times+1,4))
xi_V_B4 = np.zeros((N_times+1,4))

DV_B1 = np.zeros((N_times+1,4))
DV_B4 = np.zeros((N_times+1,4))

Dq_r_prev = 0.0

it = 0
t_0 = time.time()
for t in Times:
    
    # measurments
    grid_1 = gt.grid(syst)
    grid_1.dae2vi()
    grid_1.get_v()
    grid_1.get_i()

    phases = ['a','b','c','n']
    for iph in range(4):
        I_B1[it,iph] = float(np.abs(grid_1.I_lines[iph]))
        I_B4[it,iph] = float(np.abs(grid_1.I_lines[iph+8]))   

        
    # control
    for iph in range(3):
        I_B14 = I_B1[it,iph] - I_B4[it,iph]
        xi_V_B1[it+1,iph] = xi_V_B1[it,iph] + I_B14
        xi_V_B4[it+1,iph] = xi_V_B4[it,iph] - I_B14
        
        DV_B1[it+1,iph] = -0.0001*I_B14 -0.002*xi_V_B1[it+1,iph]
        DV_B4[it+1,iph] =  0.0001*I_B14 -0.002*xi_V_B4[it+1,iph]

    
    gt.set_voltages(syst,'B1',[V_b+DV_B1[it,0],
                               V_b+DV_B1[it,1],
                               V_b+DV_B1[it,2]],0)
                               
    gt.set_voltages(syst,'B4',[V_b+DV_B4[it,0],
                               V_b+DV_B4[it,1],
                               V_b+DV_B4[it,2]],0)
    
    events = [{'t_end':t}]
    syst.run(events)
    
#    U_grid += [U_grid_pu]
#    U += [U_meas_pu]
#    Dq_r_list += [Dq_r]

    it +=1
syst.post();    

print(time.time()-t_0)


# In[12]:


plt.close('all')
fig, axes = plt.subplots(nrows=1,ncols=2, figsize=(8, 3), frameon=False, dpi=100, squeeze=False)

for iph in range(3):
    axes[0,0].plot(Times, DV_B1[:-1,iph], label="$DV_{B1}$ (LV)")
    axes[0,1].plot(Times, DV_B4[:-1,iph], label="$DV_{B1}$ (LV)")


# In[14]:


plt.close('all')
fig, axes = plt.subplots(nrows=3,ncols=2, figsize=(10, 6), frameon=False, dpi=100, squeeze=False)

phases = ['a','b','c','n']
for iph in range(3):
    axes[iph,0].plot(Times, I_B1[:,iph], label="$DV_{B1}$ (LV)")
    axes[iph,0].plot(Times, I_B4[:,iph], label="$DV_{B1}$ (LV)")
    
    ph = phases[iph]
    
    axes[iph,1].plot(syst.T, syst.get_values(f'v_B1_{ph}_m'), label="$DV_{B1}$ (LV)")
    axes[iph,1].plot(syst.T, syst.get_values(f'v_B4_{ph}_m'), label="$DV_{B1}$ (LV)")
    
    axes[iph,0].grid()
    axes[iph,1].grid()

    axes[iph,0].plot(Times, I_B1[:,iph]*0+250, label="$DV_{B1}$ (LV)")
    axes[iph,1].plot(Times, I_B1[:,iph]*0+400/np.sqrt(3), label="$DV_{B1}$ (LV)")
    
    axes[iph,0].set_xlim([0,4])
    axes[iph,1].set_xlim([0,4])
    
    axes[iph,0].set_ylim([50,350])    
    axes[iph,1].set_ylim([400/np.sqrt(3)*0.99,400/np.sqrt(3)*1.01])    
    
    


# In[72]:


grid_1.dae2vi()
grid_1.get_v()
grid_1.get_i()
grid_1.bokeh_tools()

gt.plot_results(grid_1)


# In[73]:


syst.report_z()


# In[13]:


grid_1.V_node * grid_1.I_node


# In[7]:


310*690


# In[ ]:




