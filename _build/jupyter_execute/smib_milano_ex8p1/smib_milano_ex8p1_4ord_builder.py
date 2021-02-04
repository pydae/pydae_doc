#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sympy as sym
import numba
import pydae.build as db


# In[ ]:





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
# ### Auxiliar  equations
# 
# \begin{eqnarray}
# v_d &=& v_1\sin(\delta - \theta_1) \\
# v_q &=& v_1\cos(\delta - \theta_1) 
# \end{eqnarray}
# 
# \begin{eqnarray}
# p_e = i_d \left(v_d + R_a i_d\right) + i_q \left(v_q + R_a i_q\right) 
# \end{eqnarray}
# 
# 
# ### Differential  equations
# 
# $$
# \begin{eqnarray}
#  f_1 = \dot \delta &=&  \Omega_b \left( \omega -1 \right) \\
#  f_2 = \dot \omega &=&  \frac{1}{2H} \left( p_m - p_e - D   \left( \omega - 1 \right) \right)  \\
#  f_3 = \dot e_q &=& \frac{1}{T'_{d0}} \left( -e'_q - \left(X_d - X'_d \right) i_d + v_f^\star \right) \\
#  f_4 = \dot e_d &=&  \frac{1}{T'_{q0}} \left( -e'_d - \left(X_q - X'_q \right) i_q \right)       
# \end{eqnarray}
# $$
# 
# ### Algebraic equations
# 
# $$
# \begin{eqnarray}
# g_1 &=& v_q + R_a i_q + X_d' i_d - e_q \\
# g_2 &=& v_d + R_a i_d - X_q' i_q - e_d \\
# g_3 &=& p_t - \left(v_1 v_0 \sin \left(\theta_1 - \theta_0 \right) \right)/X_l \\
# g_4 &=& q_t + \left(v_1 v_0 \cos \left(\theta_1 - \theta_0 \right) \right)/X_l - v_1^2/X_l \\
# g_5 &=& i_d v_d + i_q v_q - p_t \\
# g_6 &=& i_d v_q - i_q v_d - q_t  
# \end{eqnarray}
# $$  
# 
# 
# 
# 
# 
# The dynamic states are:
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
# g_2\\
# g_3\\
# g_4\\
# g_5\\
# g_6
# \end{array}
# \right]
# \;\;\;\;\;\;
# $$
# 
# $$
# \mathbf x = \left[
# \begin{array}{c} 
# \delta \\ 
# \omega \\ 
#  e_q'\\ 
#  e_d' 
# \end{array} \right]
# \;\;\;\;
# \mathbf {y^{ini}} = \left[
# \begin{array}{c} 
# i_d\\
# i_q\\
# v_1\\
# \theta_1\\
# p_m\\
# v_f 
# \end{array} \right] 
# \;\;\;\;
# \mathbf {y^{run}} = \left[
# \begin{array}{c} 
# i_d\\
# i_q\\
# v_1\\
# \theta_1\\
# p_t\\
# q_t 
# \end{array} \right] 
# \;\;\;\;
# \mathbf {u^{ini}} = \left[
# \begin{array}{c} 
#  p_t\\
#  q_t 
# \end{array} \right] 
# \;\;\;\;
# \mathbf {u^{run}} = \left[
# \begin{array}{c} 
# p_m \\
# v_f
# \end{array} \right]
# $$
# 
# ### Outputs
# 
# 
# $$
# \mathbf{h} =
# \left[
# \begin{array}{c}
# p_e \\
# p_m
# \end{array}
# \right]
# \;\;\;\;\;\;
# \mathbf{z} =
# \left[
# \begin{array}{c}
# p_e\\
# p_m
# \end{array}
# \right]
# $$
#     
#     
# 
# 
#     
#     

# ## System definition 

# In[2]:


params_dict = {'X_d':1.81,'X1d':0.3,'T1d0':8.0,
               'X_q':1.76,'X1q':0.65,'T1q0':1.0,
               'R_a':0.003,'X_line': 0.05, 
               'H':3.5,'D':1.0,
               'Omega_b':2*np.pi*50,'omega_s':1.0,
               'v_0':1.0,'theta_0':0.0, 'B_shunt':0.0}


u_ini_dict = {'p_t':0.8, 'q_t':0.2}  # for the initialization problem
u_run_dict = {'p_m':0.8,'v_f':1.0}  # for the running problem (here initialization and running problem are the same)


x_list = ['delta','omega','e1q','e1d']    # [inductor current, PI integrator]
y_ini_list = ['i_d','i_q','v_1','theta_1','p_m','v_f'] # for the initialization problem
y_run_list = ['i_d','i_q','v_1','theta_1','p_t','q_t'] # for the running problem (here initialization and running problem are the same)

sys_vars = {'params_dict':params_dict,
            'u_list':u_run_dict,
            'x_list':x_list,
            'y_list':y_run_list}

exec(db.sym_gen_str())  # exec to generate the required symbolic varables and constants


# In[3]:


# Auxiliar equations
v_d = v_1*sin(delta - theta_1) 
v_q = v_1*cos(delta - theta_1) 

p_e = i_d*(v_d + R_a*i_d) + i_q*(v_q + R_a*i_q) 

# Differential equations
ddelta = Omega_b*(omega - omega_s)
domega = 1/(2*H)*(p_m - p_e - D*(omega - omega_s))
de1q = 1/T1d0*(-e1q - (X_d - X1d)*i_d + v_f)
de1d = 1/T1q0*(-e1d + (X_q - X1q)*i_q)

# Algebraic equations
g_1 = v_q + R_a*i_q + X1d*i_d - e1q
g_2 = v_d + R_a*i_d - X1q*i_q - e1d
g_3 = p_t - (v_1*v_0*sin(theta_1 - theta_0))/X_line
g_4 = q_t + (v_1*v_0*cos(theta_1 - theta_0))/X_line - v_1**2/X_line - v_1**2*B_shunt
g_5 = i_d*v_d + i_q*v_q - p_t
g_6 = i_d*v_q - i_q*v_d - q_t

# Outputs 
h_1 = p_m
h_2 = p_e

# System dictionary
sys = {'name':'smib_milano_ex8p1_4ord',
       'params_dict':params,
       'f_list':[ddelta,domega,de1q,de1d],
       'g_list':[g_1,g_2,g_3,g_4,g_5,g_6],
       'x_list':x_list,
       'y_ini_list':y_ini_list,
       'y_run_list':y_run_list,
       'u_run_dict':u_run_dict,
       'u_ini_dict':u_ini_dict,
       'h_dict':{'p_m':h_1,'p_e':h_2}}


sys = db.system(sys)
db.sys2num(sys)      # building system module
                


# In[ ]:




