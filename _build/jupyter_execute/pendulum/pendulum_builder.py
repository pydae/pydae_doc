#!/usr/bin/env python
# coding: utf-8

# \begin{equation}
# \begin{split}   \nonumber
# \mathbf {\dot x}  &   =  \mathbf{f (x,y,u) } \\
# \mathbf 0 &   =  \mathbf{g (x,y,u) }  \\
# \mathbf z &   =  \mathbf{h (x,y,u) } 
# \end{split}
# \end{equation}
# 
# - $ \mathbf x$: dynamic states 
# - $ \mathbf y$: algebraic states 
# - $ \mathbf u$: known inputs		
# - $ \mathbf f$: differential equations
# - $ \mathbf g$: algebraic equations		
# 
# 
# \begin{equation} \sf  \nonumber
# \mathbf x = \left[
# \begin{array}{c} 
# x_1 \\ 
# x_2 \\ 
# \vdots \\ 
#  x_{N_x}
# \end{array} \right] \;\;
# \mathbf y = \left[
# \begin{array}{c} 
#  y_1 \\ 
#  y_2 \\ 
# \vdots \\ 
# y_{N_y}
# \end{array} \right]
# \end{equation}
# 
# 
# 
# \begin{equation} \sf  \nonumber
# \mathbf f = \left[
# \begin{array}{c} 
#  f_1 \\ 
# f_2 \\ 
# \vdots \\ 
# f_{N_x}
# \end{array} \right] \;\;
# \mathbf g = \left[
# \begin{array}{c} 
#  g_1 \\ 
#  g_2 \\ 
# \vdots \\ 
#  g_{N_y}
# \end{array} \right]
# \end{equation}
# 
# Problema de flujo de potencia e inicialización hacia atrás:
# 	\begin{equation}
# 	\begin{split} \sf \nonumber
# 	\mathbf {\dot x}  & \sf =  \mathbf {f (x,y^{ini},u^{ini}) } \\
# 	\mathbf 0 & \sf =  \mathbf {g (x,y^{ini},u^{ini}) }  
# 	\end{split}
# 	\end{equation}
# 
# Problema de solución hacia adelante:
# 	\begin{equation}
# 	\begin{split} \sf \nonumber
# 	\mathbf {\dot x}  & \sf =  \mathbf {f (x,y^{run},u^{run}) } \\
# 	\mathbf 0 & \sf =  \mathbf {g (x,y^{run},u^{run}) }  
# 	\end{split}
# 	\end{equation}
#     
#     
# \begin{eqnarray}
# f_1 &=& \frac{dp_x}{dt} = v_x \\
# f_2 &=& \frac{dp_y}{dt} = v_y \\
# f_3 &=& \frac{dv_x}{dt} = \frac{1}{M} \left(-2 p_x \lambda + f_x - K_d v_x \right)  \\
# f_4 &=& \frac{dv_y}{dt} = \frac{1}{M} \left(-M G - 2 p_y \lambda - K_d v_y \right) \\   
# &&\\
# g_1 &=& p_x^2 + p_y^2 - L^2  \\
# g_2 &=& -\theta + \arctan\left(p_x,-p_y\right) \\
# \end{eqnarray}
# 
# The dynamic states are:
# 
# \begin{equation}
# \mathbf{f} =
# \left[
# \begin{array}{c}
# f_1\\
# f_2\\
# f_3\\
# f_4
# \end{array}
# \right]
# \end{equation}
# 
# \begin{equation} \sf  \nonumber
# \mathbf x = \left[
# \begin{array}{c} 
# p_x \\ 
# p_y \\ 
# v_x \\ 
# v_y 
# \end{array} \right] \;\;
# \mathbf {y^{ini}} = \left[
# \begin{array}{c} 
#  \lambda \\ 
#  f_x
# \end{array} \right] \;\;
# \mathbf {y^{run}} = \left[
# \begin{array}{c} 
#  \lambda \\ 
# \theta
# \end{array} \right] \;\;
# \mathbf {u^{ini}} = \left[
# \begin{array}{c} 
# \theta
# \end{array} \right] \;\;
# \mathbf {u^{run}} = \left[
# \begin{array}{c} 
# f_x
# \end{array} \right]
# \end{equation}
# 

# In[1]:


import numpy as np
import sympy as sym
import pydae.build as db


# ### Definition of variables and constants

# In[ ]:


params_dict = {'L':5.21,'G':9.81,'M':10.0,'K_d':1e-6}


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

# In[8]:


dp_x = v_x
dp_y = v_y
dv_x = (-2*p_x*lam + f_x - K_d*v_x)/M
dv_y = (-M*G - 2*p_y*lam - K_d*v_y)/M   

g_1 = p_x**2 + p_y**2 - L**2 -lam*1e-6
g_2 = -theta + sym.atan2(p_x,-p_y)


# ## Build the model

# In[9]:


sys = {'name':'pendulum',
       'params_dict':params_dict,
       'f_list':[dp_x,dp_y,dv_x,dv_y],
       'g_list':[g_1,g_2],
       'x_list':x_list,
       'y_ini_list':y_ini_list,
       'y_run_list':y_run_list,
       'u_run_dict':u_run_dict,
       'u_ini_dict':u_ini_dict,
       'h_dict':{'g_1':g_1,'PE':M*G*p_y,'KE':0.5*M*(v_x**2+v_y**2),'theta':theta}}

sys = db.system(sys)
db.sys2num(sys)


# \begin{gather*}
# a_1=b_1+c_1\\
# a_2=b_2+c_2-d_2+e_2
# \end{gather*}
# 
# \begin{align}
# a_{11}& =b_{11}&
#   a_{12}& =b_{12}\\
# a_{21}& =b_{21}&
#   a_{22}& =b_{22}+c_{22}
# \end{align}

# In[ ]:




