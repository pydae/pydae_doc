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
from pendulum import pendulum_class


# In[4]:


p = pendulum_class()


# In[5]:


M = 30.0
L = 5.21
p.initialize([{'f_x':0,'M':M,'L':L,'theta':np.deg2rad(0)}],-1)
p.report_x()
p.report_u()


# In[6]:


ssa.eval_A(p)
eig_df=ssa.damp_report(p)
eig_df


# In[7]:


freq = eig_df['Freq.']['Mode 1']
period = 1/freq
print(f'Oscillation period from small signal analysis: T = {period:0.2f} s')


# In[8]:


p.simulate([{'t_end':1, 'theta':np.deg2rad(-5)},
            {'t_end':50,'f_x':0}],'prev');


# In[9]:


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

