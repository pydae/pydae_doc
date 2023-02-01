#!/usr/bin/env python
# coding: utf-8

# ### Differential Algebraic Equations (DAE)
# 
# $$
# \begin{split}   \nonumber
# \mathbf {\dot x}  &   =  \mathbf{f (x,y,u) } \\
# \mathbf 0 &   =  \mathbf{g (x,y,u,p) }  \\
# \mathbf z &   =  \mathbf{h (x,y,u,p) } 
# \end{split}
# $$
# 
# with:
# 
# - $ \mathbf x$: dynamic states
# - $ \mathbf y$: algebraic states
# - $ \mathbf u$: known inputs	
# - $ \mathbf p$: parameters		
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
# Backward solution:
# 
# $$
# \begin{split} 
# \mathbf {\dot x}  &  =  \mathbf {f (x,y^{ini},u^{ini},p) } \\
# \mathbf 0 & =  \mathbf {g (x,y^{ini},u^{ini},p) }  
# \end{split}
# $$
# 
# Foreward solution:
# 
# $$
# \begin{split} 
# \mathbf {\dot x}  &  =  \mathbf {f (x,y^{run},u^{run},p) } \\
# \mathbf 0 & =  \mathbf {g (x,y^{run},u^{run},p) }  
# \end{split}
# $$
# 
# 

# 
