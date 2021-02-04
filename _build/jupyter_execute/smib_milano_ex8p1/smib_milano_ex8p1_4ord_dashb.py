#!/usr/bin/env python
# coding: utf-8

# # SMIB system as in Milano's book example 8.1 (4ord order added)

# In[1]:


get_ipython().run_line_magic('matplotlib', 'widget')


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sopt
import ipywidgets


# ## Import system module

# In[3]:


from smib_milano_ex8p1_4ord import smib_milano_ex8p1_4ord_class,daesolver


# ## Instantiate system

# In[4]:


syst = smib_milano_ex8p1_4ord_class()
syst.t_end = 15.0
syst.Dt = 0.05
syst.decimation =1
syst.update()


# ## Solution function

# In[5]:


def test(syst,p_m,v_f,use_numba=True):
    # simulation parameters
    syst.struct[0].imax = 50    # maximum number of iterations
    syst.struct[0].itol = 1e-8  # relative tolerance to stop iteration
    syst.struct[0].solvern = 5  # 5 = DAE trapezoidal solver (fixed step)
    syst.struct[0].it = 0       # set time step to zero
    syst.struct[0].it_store = 0 # set storage to zero
    syst.struct[0].t = 0.0      # set time to zero
    
    syst.struct[0].D = 0
    
    syst.struct[0].p_m = 0
    syst.struct[0].v_f = 0
    syst.struct[0].T1d0 = 4

    syst.struct.P_t = 0
    syst.struct.Q_t = 0    
    
    # compute initial conditions using x and y_ini 
    xy0 = np.ones(syst.N_x+syst.N_y)
    xy = sopt.fsolve(syst.ini_problem,xy0 )

    # from ini system to run system
    syst.struct[0].p_m = xy[syst.xy_ini_list.index('p_m')]
    syst.struct[0].v_f = xy[syst.xy_ini_list.index('v_f')]

    syst.struct[0].x[:,0] = xy[0:syst.N_x]
    syst.struct[0].y[:,0] = xy[syst.N_x:]

    syst.struct[0].y[syst.y_list.index('P_t'),0] = syst.struct.P_t
    syst.struct[0].y[syst.y_list.index('Q_t'),0] = syst.struct.Q_t
    
    # solve system
    syst.struct.t_end = 1.0  
    daesolver(syst.struct)    # run until 1 s
    syst.struct[0].p_m = p_m  # apply step in mechanical power p_m
    syst.struct[0].v_f = v_f  # apply step in mechanical power p_m
    syst.struct.t_end = 20.0  
    daesolver(syst.struct)    # run until 10 s
        
    T = syst.struct[0]['T'][:syst.struct[0].it_store]
    X = syst.struct[0]['X'][:syst.struct[0].it_store,:]
    Y = syst.struct[0]['Y'][:syst.struct[0].it_store,:]

    return T,X,Y


# ## Run test

# In[6]:


T,X,Y =test(syst,0.0,1.0,use_numba=True)


# In[ ]:





# ### Results extraction

# In[16]:


plt.ioff()
plt.clf()

T,X,Y =test(syst,0.0,1.0,use_numba=True)
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(5, 4), frameon=False)
fig.canvas.toolbar_visible = False

line_delta = axes[0,0].plot(T, X[:,syst.x_list.index('delta')], label='$\sf \delta$')
line_omega = axes[1,0].plot(T, X[:,syst.x_list.index('omega')], label='$\sf \omega$')
line_v_1 = axes[0,1].plot(T, Y[:,syst.y_list.index('v_1')], label='$\sf v_1$')
#line_theta_1 = axes[0,1].plot(T, Y[:,syst.y_list.index('theta_1')], label='$\sf \\theta_1$')
line_p_t = axes[1,1].plot(T, Y[:,syst.y_list.index('P_t')], label='$\sf P_t$')
line_q_t = axes[1,1].plot(T, Y[:,syst.y_list.index('Q_t')], label='$\sf Q_t$')

x_0 = X[0,:]
y_0 = Y[0,:]

y_labels = ['$\delta$','$\omega$','$P_t$']

axes[0,0].set_ylim((-1,2))
axes[1,0].set_ylim((0.95,1.05))
axes[0,1].set_ylim((0.8,1.2))
axes[1,1].set_ylim((-0.5,1.5))

axes[0,0].grid(True)
axes[1,0].grid(True)
axes[0,1].grid(True)
axes[1,1].grid(True)
axes[0,0].legend(loc='best')
axes[1,0].legend(loc='best')
axes[0,1].legend(loc='best')
axes[1,1].legend(loc='best')

axes[1,0].set_xlabel('Time (s)')  
axes[1,1].set_xlabel('Time (s)') 

fig.tight_layout()
#axes[0].set_title('Par en función de la velocidad')
#axes[1].set_title('Corriente en función de la velocidad')


sld_p_m = ipywidgets.FloatSlider(orientation='horizontal',description = u"p_m:", 
                                value=0.0, min=0.0,max= 1.2, 
                                step=.1)


sld_v_f = ipywidgets.FloatSlider(orientation='horizontal',description = u"v_f:", 
                                value=syst.struct.v_f, min=0.0,max= 3.0, 
                                step=.1)

prog_c = ipywidgets.IntProgress(
    value=100,
    min=0,
    max=120,
    step=1,
    description='Carga:',
    bar_style='', # 'success', 'info', 'warning', 'danger' or ''
    orientation='horizontal' 
)


# ### Results plots

# In[17]:




def update(change):
   
   p_m = sld_p_m.value
   v_f = sld_v_f.value
   
   try:
       T,X,Y = test(syst,p_m,v_f,use_numba=True)

       line_delta[0].set_data(T, X[:,syst.x_list.index('delta')])
       line_omega[0].set_data(T, X[:,syst.x_list.index('omega')])
       line_v_1[0].set_data(T, Y[:,syst.y_list.index('v_1')])
       #line_theta_1[0].set_data(T, Y[:,syst.y_list.index('theta_1')])
       line_p_t[0].set_data(T, Y[:,syst.y_list.index('P_t')])
       line_q_t[0].set_data(T, Y[:,syst.y_list.index('Q_t')])

       c = np.abs(Y[-1,syst.y_list.index('i_d')]+1j*Y[-1,syst.y_list.index('i_q')])

       prog_c.bar_style = 'success'
       if c>0.9:
           prog_c.bar_style = 'warning'
       if c>1.0:
           prog_c.bar_style = 'danger'
       prog_c.value = 100*c
       
       fig.canvas.draw_idle()
   except:
         print("An exception occurred")        
  
   
sld_p_m.observe(update, names='value')
sld_v_f.observe(update, names='value')

layout_row1 = ipywidgets.HBox([fig.canvas])
layout_row2 = ipywidgets.HBox([sld_p_m,sld_v_f,prog_c])

layout = ipywidgets.VBox([layout_row1,layout_row2])
layout


# In[ ]:





# In[ ]:




