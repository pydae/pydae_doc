#!/usr/bin/env python
# coding: utf-8

# # Pendulum

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter  
from IPython.core.display import HTML


# ## Import system module

# In[3]:


from pydae import ssa
import pendulum


# In[4]:


p = pendulum.model()


# In[5]:


M = 30.0
L = 5.21
p.ini({'f_x':0,'M':M,'L':L,'theta':np.deg2rad(0)},-1)
p.report_x()
p.report_u()


# In[6]:


ssa.A_eval(p)
eig_df=ssa.damp_report(p)
eig_df.round(3)


# In[7]:


freq = eig_df['Freq.']['Mode 1']
period = 1/freq
print(f'Oscillation period from small signal analysis: T = {period:0.2f} s')


# In[8]:


p.ini({'f_x':0,'M':M,'L':L,'theta':np.deg2rad(-5)},-1)
p.run( 1.0, {})
p.run(50.0, {'f_x':0}) 
p.post();


# In[9]:


time = p.Time
theta = np.rad2deg(p.get_values('theta'))

idx_1 = np.where(theta==np.max(theta[(time>7)&(time<11)]))[0][0]
idx_2 = np.where(theta==np.max(theta[(time>10)&(time<14)]))[0][0]
t_1 = time[idx_1]
t_2 = time[idx_2]

period_sim = t_2 - t_1

print(f'Oscillation period from simulation: T = {period_sim:0.2f} s')

plt.close('all')
fig, axes = plt.subplots(nrows=1,ncols=1, figsize=(5, 3), dpi=100)

axes.plot(p.Time, np.rad2deg(p.get_values('theta')), label=f'$\theta$')
axes.grid()
axes.set_ylabel('$\\theta (º)$')


# In[10]:


times = p.Time
t_end = times[-1] 
pos = p.get_values('theta')
keyTimes = ""
keyPoints = ""
for it in range(len(times)):
    keyTimes  += f'{times[it,0]/t_end:0.4f};'
    keyPoints += f'{(-pos[it]+0.5):0.4f};' 
keyTimes  = keyTimes[:-1].replace("'",'"')    
keyPoints = keyPoints[:-1].replace("'",'"') 

#keyPoints = "0;0.2;0.4;0.6;0.8;1.0"
#keyTimes  = "0;0.3;0.5;0.6;0.8;1.0"

HTML(f'''<svg viewBox="0 0 500 500" xmlns="http://www.w3.org/2000/svg"
    xmlns:xlink="http://www.w3.org/1999/xlink">
  <path d="m 110.36905,102.05357 a 39.309525,39.309525 0 0 1 -19.654764,34.04305 39.309525,39.309525 0 0 1 -39.309525,0 A 39.309525,39.309525 0 0 1 31.75,102.05357"
      stroke="lightgrey" stroke-width="2" fill="none" id="motionPath"/>
  <circle r="5" fill="red">
    <animateMotion dur="{t_end}" repeatCount="indefinite"
        keyPoints={keyPoints} 
        keyTimes= {keyTimes} calcMode="linear">
      <mpath xlink:href="#motionPath"/>
    </animateMotion>
  </circle>
</svg>''')


# In[32]:


times = p.T
t_end = times[-1,0] 
pos = np.rad2deg(p.get_values('theta'))
keyTimes = ""
keyPoints = ""
for it in range(len(times)):
    keyTimes  += f'{times[it,0]/t_end:0.4f};'
    keyPoints += f'{(-pos[it]+0.5):0.4f},71.4375,103.18749;' 
keyTimes  = keyTimes[:-1].replace("'",'"')    
keyPoints = keyPoints[:-1].replace("'",'"') 


HTML(f'''
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   id="svg8"
   version="1.1"
   viewBox="0 0 104.56544 58.339157"
   height="58.339157mm"
   width="104.56544mm">
  <defs
     id="defs2" />
  <metadata
     id="metadata5">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     id="layer1"
     transform="translate(-19.155907,-101.33694)">
    <path
       style="fill:none;stroke:#555555;stroke-width:0.5;stroke-miterlimit:4;stroke-dasharray:1, 1;stroke-dashoffset:0;stop-color:#000000"
       id="path833"
       d="m 123.47135,103.1875 a 52.033852,52.033405 0 0 1 -26.016925,45.06225 52.033852,52.033405 0 0 1 -52.033852,0 52.033852,52.033405 0 0 1 -26.016925,-45.06225" />
    <g transform="rotate(0,71.4375,103.18749)">
       id="g890">
      <rect
         style="fill:#337ab7;stroke:none;stroke-width:0.499999;stroke-miterlimit:4;stroke-dasharray:0.999999, 0.999999;stroke-dashoffset:0;stop-color:#000000"
         id="rect885"
         width="1.5874945"
         height="50.270844"
         x="70.643753"
         y="101.86458" />
      <circle
         style="fill:#d9534f;stroke:#d9534f;stroke-width:0.132292;stop-color:#000000"
         id="path868"
         cx="71.4375"
         cy="155.25372"
         r="4.3562231" />
      <circle
         style="fill:#337ab7;stroke:none;stroke-width:0.132292;stop-color:#000000"
         id="path868-1"
         cx="71.4375"
         cy="103.18751"
         r="1.8505714" />
    <animateTransform attributeType="xml" attributeName="transform" type="rotate"
    values={keyPoints}
        keyTimes={keyTimes} 
        dur="{t_end}s" repeatCount="indefinite" />
    </g>
  </g>
</svg>
''')


# ### Animation

# In[10]:


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



