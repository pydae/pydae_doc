#!/usr/bin/env python
# coding: utf-8

# ### Differential Algebraic Equations (DAE)
# 
# $$
# \begin{split}   \nonumber
# \mathbf {\dot x}  &   =  \mathbf{f (x,y,u) } \\
# \mathbf 0 &   =  \mathbf{g (x,y,u) }  \\
# \mathbf z &   =  \mathbf{h (x,y,u) } 
# \end{split}
# $$
# 
# with:
# 
# - $ \mathbf x$: dynamic states
# - $ \mathbf y$: algebraic states
# - $ \mathbf u$: known inputs		
# - $ \mathbf f$: differential equations
# - $ \mathbf g$: algebraic equations	
# - $ \mathbf z$: outputs
# - $ \mathbf h$: outputs functions
# 
# 
# $$ 
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
# $$
# 
# 
# 
# $$
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
# $$
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
# $$   
# \begin{eqnarray}
# f_1 &=& \frac{dp_x}{dt} = v_x \\
# f_2 &=& \frac{dp_y}{dt} = v_y \\
# f_3 &=& \frac{dv_x}{dt} = \frac{1}{M} \left(-2 p_x \lambda + f_x - K_d v_x \right)  \\
# f_4 &=& \frac{dv_y}{dt} = \frac{1}{M} \left(-M G - 2 p_y \lambda - K_d v_y \right) \\   
# &&\\
# g_1 &=& p_x^2 + p_y^2 - L^2  \\
# g_2 &=& -\theta + \arctan\left(p_x,-p_y\right) \\
# \end{eqnarray}
# $$
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
# $$
# 
# $$
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
# $$
# 
