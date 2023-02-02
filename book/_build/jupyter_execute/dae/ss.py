#!/usr/bin/env python
# coding: utf-8

# ## Steady state solution
# 
# The steady state of a dynamical system is obtained when the derivatives are zero ($\mathbf {\dot x = \mathbf 0}$). This can be extended to DAE system as follows:
# 
# $$
# \begin{split} 
# \mathbf 0 & =  \mathbf {f (x,y,u,p) } \\
# \mathbf 0 & =  \mathbf {g (x,y,u,p) }  
# \end{split}
# $$
# 
# In pydae the steady state solution is obtained for the initialization problem:
# 
# 
# $$
# \begin{split} 
# \mathbf 0 & =  \mathbf {f (x,y^{ini},u^{ini},p) } \\
# \mathbf 0 & =  \mathbf {g (x,y^{ini},u^{ini},p) }  
# \end{split}
# $$

# Therefore, to find the steady state one must solve the previous system for $\mathbf x$ and $\mathbf{ y^{ini}}$ or for the vector containing both called $\boldsymbol \xi$:
# 
# $$
# \boldsymbol \lambda \left( \boldsymbol \xi \right) =  \left[
# \begin{array}{c}
# \mathbf {f (\boldsymbol \xi)}	\\ 
# \mathbf {g (\boldsymbol \xi)}
# \end{array} \right] \;\;\;\;\;\;
# \boldsymbol \xi  =  \left[
# \begin{array}{c}
# \mathbf {x}	\\ 
# \mathbf {{y}^{ini}}
# \end{array} \right]
# $$
# 
# for given inputs $\mathbf u^{ini}$ and parameters $\mathbf p$.
# 

# Thus, to find the equilibrium point, the following system of $ N_x + N_y$ equations must be solved:
# 
# $$
# \mathbf 0 = \boldsymbol \lambda \left( \boldsymbol \xi \right) 
# $$

# In pydae this system is can be solved with Newton Raphson method:
# 
# $$ \Delta \boldsymbol {\xi} =-\left( \left. \mathbf{J_{ini}} \right|_{\boldsymbol {\xi} ^i} \right)^{-1}  \boldsymbol \lambda \left( \boldsymbol \xi^i \right)  
# $$
# 
# To find a new value of $\boldsymbol {\xi}$:
# 
# $$ \boldsymbol \xi^{i+1} =\boldsymbol \xi^i +  \Delta \boldsymbol {\xi}   $$
# 
# where the jacobian is defined as follows:
# 
# $$   \mathbf {J_{ini} }= \left[ 
# \begin{array}{cc}
# \mathbf{F_x^{ini}} & \mathbf{F_y^{ini}}	\\ 
# \mathbf{G_x^{ini}} & \mathbf{G_y^{ini}}	 
# \end{array}			\right]
# $$

# In order to improve the solution speed, the jacobian $\mathbf{ J_{ini}}$ is composed by three matrices as follows:
# $$ \mathbf {J_{ini} =  J_{ini}^{num} +  J_{ini}^{up} + J_{ini}^{xy} } $$  
# where:
# 
# * $\mathbf {J_{ini}^{num}}$ has nonzero elements that are constants and is only evaluated when the system is build.
# * $\mathbf {J_{ini}^{up}}$ has nonzero elements that are only function of the inputs $\mathbf {u^{ini}}$ or the parameters   $\mathbf p$. These elements are only evaluated when inputs and parameters are changed.
# * $\mathbf {J_{ini}^{xy}}$ has nonzero elements that contains dynamics or algebraic states $\mathbf x$ and $\mathbf {y^{ini}}$.
# 

# In[ ]:




