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
