#!/usr/bin/env python
# coding: utf-8

# # Pendulum DAE system

# In[1]:


from IPython.core.display import HTML
from IPython.lib.display import YouTubeVideo
YouTubeVideo('4a0FbQdH3dY')


# - 9.43 he sets the bom in 5deg (How much force in x did he applied?)
# - 11.35 he sets the bom in 10deg (How much force in x did he applied?)
# - 16:20 change in mass
# - 24:00 without little push  
# - 24:42 with little push glass broken
# - 26:50 dagerous experiment

# In[2]:


HTML('../svg/pendulum_doc.svg')


# ## Formulation
# 
# Backward solution:
# 
# $$
# \begin{split} 
# \mathbf {\dot x}  &  =  \mathbf {f (x,y^{ini},u^{ini}) } \\
# \mathbf 0 & =  \mathbf {g (x,y^{ini},u^{ini}) }  
# \end{split}
# $$
# 
# Foreward solution:
# 
# $$
# \begin{split} 
# \mathbf {\dot x}  &  =  \mathbf {f (x,y^{run},u^{run}) } \\
# \mathbf 0 & =  \mathbf {g (x,y^{run},u^{run}) }  
# \end{split}
# $$
#  
# ### Differential equations
# 
# $$   
# \begin{eqnarray}
# f_1 &=& \frac{dp_x}{dt} = v_x \\
# f_2 &=& \frac{dp_y}{dt} = v_y \\
# f_3 &=& \frac{dv_x}{dt} = \frac{1}{M} \left(-2 p_x \lambda + f_x - K_d v_x \right)  \\
# f_4 &=& \frac{dv_y}{dt} = \frac{1}{M} \left(-M G - 2 p_y \lambda - K_d v_y \right)
# \end{eqnarray}
# $$
# 
# ### Algebraic equations
# 
# $$
# \begin{eqnarray}
# g_1 &=& p_x^2 + p_y^2 - L^2  \\
# g_2 &=& -\theta + \arctan\left(p_x,-p_y\right) 
# \end{eqnarray}
# $$
# 
# ### DAE equations in vector form
# 
# $$
# \mathbf{f} =
# \left[
# \begin{array}{c}
# f_1\\
# f_2\\
# f_3\\
# f_4
# \end{array}
# \right]
# \;\;\;\;\;\;
# \mathbf{g} =
# \left[
# \begin{array}{c}
# g_1\\
# g_2
# \end{array}
# \right]
# \;\;\;\;\;\;
# $$
# 
# ### Dynamic and algebraic states and inputs
# 
# $$
# \mathbf x = \left[
# \begin{array}{c} 
# p_x \\ 
# p_y \\ 
# v_x \\ 
# v_y 
# \end{array} \right]
# \;\;\;\;
# \mathbf {y^{ini}} = \left[
# \begin{array}{c} 
#  \lambda \\ 
#  f_x
# \end{array} \right] 
# \;\;\;\;
# \mathbf {y^{run}} = \left[
# \begin{array}{c} 
#  \lambda \\ 
# \theta
# \end{array} \right] 
# \;\;\;\;
# \mathbf {u^{ini}} = \left[
# \begin{array}{c} 
# \theta
# \end{array} \right] 
# \;\;\;\;
# \mathbf {u^{run}} = \left[
# \begin{array}{c} 
# f_x
# \end{array} \right]
# $$
# 
# ### Outputs
# 
# $$
# \begin{split} 
# \mathbf z & =  \mathbf {h (x,y^{run},u^{run}) }  
# \end{split}
# $$
# 
# We may be interested in the potential energy:
# 
# $$
# E_p = M G \left(p_y + L\right)
# $$
# 
# And the kinetical energy:
# 
# $$
# E_k = \frac{1}{2} M \left(v_x^2 + v_y^2 \right)
# $$
# 
# $$
# \mathbf{h} =
# \left[
# \begin{array}{c}
# M G p_y\\
# \frac{1}{2} M \left(v_x^2 + v_y^2 \right)\\
# f_x
# \end{array}
# \right]
# \;\;\;\;\;\;
# \mathbf{z} =
# \left[
# \begin{array}{c}
# E_p\\
# E_k\\
# f_x
# \end{array}
# \right]
# $$

# ## pydae model building

# In[3]:


import numpy as np
import sympy as sym
import pydae.build as db


# ### Definition of variables and constants

# In[4]:


params_dict = {'L':5.21,'G':9.81,'M':10.0,'K_d':1e-3}  # parameters with default values


u_ini_dict = {'theta':np.deg2rad(5.0)}  # input for the initialization problem
u_run_dict = {'f_x':0}                  # input for the running problem, its value is updated 


x_list = ['p_x','p_y','v_x','v_y']  # dynamic states
y_ini_list = ['lam','f_x']          # algebraic states for the initialization problem
y_run_list = ['lam','theta']        # algebraic for the running problem

sys_vars = {'params':params_dict,
            'u_list':u_run_dict,
            'x_list':x_list,
            'y_list':y_run_list}

exec(db.sym_gen_str())  # exec to generate the required symbolic varables and constants


# ### System formulation

# In[5]:


dp_x = v_x
dp_y = v_y
dv_x = (-2*p_x*lam + f_x - K_d*v_x)/M
dv_y = (-M*G - 2*p_y*lam - K_d*v_y)/M   

g_1 = p_x**2 + p_y**2 - L**2 -lam*1e-6
g_2 = -theta + sym.atan2(p_x,-p_y)


# ### Build the model

# In[6]:


sys = {'name':'pendulum',
       'params_dict':params_dict,
       'f_list':[dp_x,dp_y,dv_x,dv_y],
       'g_list':[g_1,g_2],
       'x_list':x_list,
       'y_ini_list':y_ini_list,
       'y_run_list':y_run_list,
       'u_run_dict':u_run_dict,
       'u_ini_dict':u_ini_dict,
       'h_dict':{'E_p':M*G*(p_y+L),'E_k':0.5*M*(v_x**2+v_y**2),'f_x':f_x}}

sys = db.system(sys)
db.sys2num(sys)


# In[ ]:





# In[7]:


HTML('''
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 200" width="400" height="200" version="1.1">

  <g>
    <!--control starts with radius 100 but shrinks to radius 0 when clicked.  It goes back to 100 when the reverseanimation runs-->
    <circle cx="100" cy="100" r="100" fill="#b2d4e5">
      <animate id="startAnimation" dur="1.0s" attributeName="r" values="100; 0" fill="freeze" begin="click" />
      <animate dur="1.0s" attributeName="r" values="0; 100" fill="freeze" begin="reverseAnimation.begin" />
    </circle>
    <text y="240" x="50" font-family="Verdana" text-align="center" font-size="30" textLength="100">Start</text>
  </g>

  <g transform="translate(300,0)">
    <circle cx="100" cy="100" r="100" fill="#b2d4e5">
    </circle>

    <!-- vertical line -->
    <path d="M100 50 l0 100" stroke="white" stroke-width="20" stroke-linecap="round">
      <animate dur="0.5s" begin="startAnimation.begin" attributeName="d" values="M100 50 l0 100; M100 100 l0 0; M100 150 l50 -100" fill="freeze" />
      <animate dur="0.5s" begin="reverseAnimation.begin" attributeName="d" values="M100 150 l50 -100; M100 100 l0 0; M100 50 l0 100" fill="freeze" />
    </path>

    <!-- horizontal line -->
    <path d="M50 100 l100 0" stroke="white" stroke-width="20" stroke-linecap="round">
      <animate dur="0.5s" begin="startAnimation.begin" attributeName="d" values="M50 100 l100 0; M100 100 l0 0; M50 100 l50 50" fill="freeze" />
      <animate dur="0.5s" begin="reverseAnimation.begin" attributeName="d" values=" M50 100 l50 50; M100 100 l0 0; M50 100 l100 0" fill="freeze" />
    </path>
  </g>

  <g transform="translate(600,0)">
    <!--control starts with radius 0 but grows to radius 100 when startanimation runs.  It goes back to 0 clicked-->
    <circle cx="100" cy="100" r="0" fill="#b2d4e5">
      <animate dur="1s" attributeName="r" values="0; 100" fill="freeze" begin="startAnimation.begin" />
      <animate id="reverseAnimation" dur="1s" attributeName="r" values="100; 0" fill="freeze" begin="click" />
    </circle>
    <text y="240" x="25" font-family="Verdana" text-align="center" font-size="30" textLength="150">Reverse</text>
  </g>
</svg>
''')


# In[ ]:




