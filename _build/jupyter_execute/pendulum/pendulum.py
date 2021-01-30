#!/usr/bin/env python
# coding: utf-8

# # Pendulum

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter  


# ## Import system module

# In[3]:


from pydae import ssa
from pendulum import pendulum_class


# In[152]:


p = pendulum_class()


# In[153]:


M = 30.0
L = 5.21
p.initialize([{'f_x':0,'M':M,'L':L,'theta':np.deg2rad(0)}],-1)
p.report_x()
p.report_u()


# In[154]:


ssa.eval_A(p)
eig_df=ssa.damp_report(p)
eig_df


# In[155]:


freq = eig_df['Freq.']['Mode 1']
period = 1/freq
print(f'Oscillation period from small signal analysis: T = {period:0.2f} s')


# In[156]:


p.simulate([{'t_end':1, 'theta':np.deg2rad(-5)},
            {'t_end':50,'f_x':0}],'prev');


# In[157]:


time = p.T[:,0]
theta = np.rad2deg(p.get_values('theta'))

idx_1 = np.where(theta==np.max(theta[(time>7)&(time<11)]))[0][0]
idx_2 = np.where(theta==np.max(theta[(time>10)&(time<14)]))[0][0]
t_1 = time[idx_1]
t_2 = time[idx_2]

period_sim = t_2 - t_1

print(f'Oscillation period from simulation: T = {period_sim:0.2f} s')

plt.close('all')
fig, axes = plt.subplots(nrows=1,ncols=1, figsize=(5, 3), dpi=100)

axes.plot(p.T, np.rad2deg(p.get_values('theta')), label=f'$\theta$')
axes.grid()
axes.set_ylabel('$\\theta (º)$')


# In[166]:


(85**(1/3))*5


# ### Animation

# In[179]:


M_list = [10,10,85]
theta_ini_list = [5,10,10]

for M,theta_ini in zip(M_list,theta_ini_list):
    
    p.simulate([{'t_end':1, 'theta':np.deg2rad(-theta_ini),'M':M},
                {'t_end':50,'f_x':0}],'prev');

    fig, ax = plt.subplots() 
    line, = ax.plot([], [], 'b', lw=2) 
    pivot, = ax.plot(0, 0, 'bo', ms=7) 
    
    bob, = ax.plot([], [], 'ro', ms=(M**(1/3))*5)  

    txt = ax.text(1,0.5,f't = 0 s')
    txt_theta = ax.text(1,0.0,f'$\\theta = 0 º')
    def init():  
        ax.set_xlim(-1.5,1.5)  
        ax.set_ylim(-1.2*L,1) 


    def update(t): 
        pos_x = np.interp(t,p.T[:,0],p.get_values('pos_x'))  
        pos_y = np.interp(t,p.T[:,0],p.get_values('pos_y')) 
        theta = np.interp(t,p.T[:,0],p.get_values('theta')) 
        bob.set_data(pos_x, pos_y)
        txt.set_text(f't = {t:5.1f} s')
        txt_theta.set_text(f'$\\theta$ = {np.rad2deg(theta):5.1f} º')
        x = [0,pos_x]  
        y = [0,pos_y] 
        line.set_data(x, y)    

    Dt_ani = 1/25    
    ani = FuncAnimation(fig, update,frames=np.arange(0,50,Dt_ani),init_func=init,interval=Dt_ani, blit=True);
    writer = PillowWriter(fps=np.ceil(1/Dt_ani))  
    ani.save(f"pendulum_{M:0.0f}kg_{theta_ini:0.0f}deg.png", writer=writer)  


# In[171]:


from IPython.display import HTML


# In[176]:


ani = FuncAnimation(fig, update,frames=np.arange(0,50,Dt_ani),init_func=init,interval=Dt_ani);
    writer = PillowWriter(fps=np.ceil(1/Dt_ani))
HTML(ani.to_html5_video())


# In[131]:


import numpy as np
import matplotlib.pyplot as plt


r = np.arange(0, 2, 0.01)
theta = 2 * np.pi * r

ax = plt.subplot(projection='polar')
ax.plot(theta, r)
ax.set_rmax(2)
ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()


# In[31]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    xdata.append(frame)
    ydata.append(np.sin(frame))
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128),
                    init_func=init, blit=True)
plt.show()


# 

# In[178]:


import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

t = np.linspace(0,2*np.pi)
x = np.sin(t)

fig, ax = plt.subplots()
ax.axis([0,2*np.pi,-1,1])
l, = ax.plot([],[])

def animate(i):
    l.set_data(t[:i], x[:i])

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=len(t))

from IPython.display import HTML
HTML(ani.to_html5())


# In[ ]:



