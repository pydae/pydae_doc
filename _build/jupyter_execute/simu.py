#!/usr/bin/env python
# coding: utf-8

# # Ejemplo de simulación

# In[1]:



1.3*250


# In[2]:


1.45*335


# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from pydae.systems.ier.proyecto import proyecto_class
import pydae.grid_tools as gt
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")


# Primero se instancia el simulador:

# In[4]:


pr = proyecto_class()


# El sistema tiene entradas, estados y salidas. Siendo los más importantes lo que se muestran en la siguiente figura:
# 
# ![title](./svg/proyecto_4cig_pydae.svg)

# Los parametros de las líneas y transformadores deben ser adaptados, para ello se puede utilizar la herramienta `change_line` del modulo `pydae.grid_tools`:

# In[5]:


gt.change_line(pr,'GR4','GR3',R_km=8.0,X_km=4.0,km=0.5)
gt.change_line(pr,'GR2','GR1',R_km=8.0,X_km=4.0,km=0.5)
gt.change_line(pr,'GR3','PMV',R_km=4.0,X_km=2.0,km=0.5)
gt.change_line(pr,'GR1','PMV',R_km=4.0,X_km=2.0,km=0.5)
gt.change_line(pr,'PMV','POI',R_pu=0.01,X_pu=0.05,S_mva=52.0)
gt.change_line(pr,'POI','GRI',R_km=0.168495,X_km=0.04,km=26)


# Una vez que el sistema tiene los parametros deseados se lo puede inicializar.

# In[6]:


P_GR = 5e6
Q_GR = 5e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR}], 1) 


# Para comprobar que la inicialización fue exitosa se pueden consultar los valores de las variables algebráicas:

# In[7]:


pr.report_y()


# `pydae.grid_tools` también aporta métodos para extraer información del sistema resuelto. Por ejemplo se puede obtener la corriente o la potencia compleja por una línea:

# In[8]:


gt.get_line_i(pr,'POI','GRI')


# In[9]:


gt.get_line_s(pr,'POI','GRI') 


# Para poder simular la planta y un controlador se pueden emplear los siguientes 3 pasos:
# 
# 1. Inicialización (initialize)
# 2. Correr hasta un determinado tiempo (run)
# 3. Postprocesar (post)
# 
# El paso 2 se puede ejecutar repetitivamante para ir cambiando las entradas del sistema a lo largo del tiempo. A continuación se simula hasta ts= 1 s, se realiza un cambio en la potencia reactiva del STATCOM (`Q_PMV`) y se simula hasta t = 10 s:

# In[10]:


pr.initialize([{'Q_PMV':0.0}],1)
print(f"V_POI antes del cambio: {pr.get_value('V_POI')}")
pr.run([{'t_end':1.0}])
pr.set_value('Q_PMV',10e6)
pr.run([{'t_end':10.0}])
print(f"V_POI después del cambio: {pr.get_value('V_POI')}")
pr.post();


# Al ser una simulación en el dominio del tiempo luego de terminar la simulación con `pr_post()`, se puede graficar la evolución de la tensión y de la potencia reactiva:

# In[11]:


fig0, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')
axes[1].plot(pr.T,pr.get_values('Q_PMV')/1e6,label='$Q_{PMV}$')
axes[0].set_ylabel('Voltage (p.u.)')
axes[1].set_ylabel('Reactive Power (Mvar)')
axes[1].set_xlabel('Time (s)')


# ## Control de tensión del POI
# 
# Una función importante de un contralador PPC es el de controlar la tensión en una determinada barra, como puede ser la denominada POI en el diagrama de arriba. Con este objetivo el PPC medirá una tensión `V_POI` y enviara la referencia de reactiva al STATCOM (`Q_PMV`) y a los generadores (`Q_GR1`, `Q_GR2`, `Q_GR3` y `Q_GR4`).
# 
# A continuación se presentan dos casos:
# 
# 1. Control de tensión utilizando solamente el STATCOM.
# 2. Control de tensión coordinado STATCOM+Generadores
# 
# ### Control de tensión utilizando sólo el STATCOM

# In[12]:


# Inicialización del sistema
# ==========================

P_GR = 5e6
Q_GR = 0e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR,
                'P_PMV': 0.0, 'Q_PMV': 0.0}], 1) 

K_g = 1/250e6 # gain of the plant

T_s = 0.1
times = np.arange(0,10,T_s)
N_t = len(times)
N_x = 2

# diseño del controlador
K_p = -20.0e6
K_i = 1/(T_s*K_g)

# condiciones iniciales
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
V_POI_ref = pr.get_value('V_POI')

# vector para guardar estados, etc.
out = np.zeros((N_t,N_x))
it = 0


# Bucle de simulación
# ===================

for t in times:
    out[it,0] = xi
    out[it,1] = Q_POI_ref
    
    if t>1.0:
        V_POI_ref = 1.05
        
    # Plant
    pr.set_value('Q_PMV',Q_POI_ref)
    pr.run([{'t_end':t,'Q_PMV':Q_POI_ref}])
    V_POI = pr.get_value('V_POI')
    
    # Controller
    error = V_POI_ref - V_POI
    xi = xi + T_s*error
    Q_POI_ref = K_p*error + K_i*xi
    
    it += 1 

# Post-procesado
# ==============
pr.post()

fig1, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].step(times,out[:,1]/1e6,label='$Q_{POI}^\star$')
axes[0].step(times,out[:,1]/1e6,label='$Q_{PMV}^\star$')

axes[1].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')

axes[0].set_ylabel('Reactive Power (Mvar)')
axes[0].set_ylim([-5,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,10])

ax.set_xlabel('Time (s)')
fig1.align_ylabels()


# In[13]:


# Inicialización del sistema
# ==========================

P_GR = 5e6
Q_GR = 0e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR,
                'P_PMV': 0.0, 'Q_PMV': 0.0}], 1) 

K_g = 1/250e6 # gain of the plant

T_s = 0.1
times = np.arange(0,10,T_s)
N_t = len(times)
N_x = 2

# diseño del controlador
K_p = -20.0e6
K_i = 1/(T_s*K_g)

# condiciones iniciales
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
V_POI_ref = pr.get_value('V_POI')
Q_STATCOM_max = 10.0e6
# vector para guardar estados, etc.
out = np.zeros((N_t,5))
it = 0
Q_STATCOM_ref = 0
Q_GR_ref = 0
# Bucle de simulación
# ===================

for t in times:
    out[it,0] = xi
    out[it,1] = Q_POI_ref
    out[it,2] = Q_STATCOM_ref
    
    if t>1.0:
        V_POI_ref = 1.05
        
    # Plant
    Q_GR1_sat = Q_GR1_ref 
    Q_GR1_max = np.sqrt(S_GR1**2 - P_GR1**2)
    if Q_GR1_sat > Q_GR1_max: Q_GR1_sat  = Q_GR1_max

    pr.run([{'t_end':t,
             'Q_PMV':Q_STATCOM_ref,
             'Q_GR1':Q_GR1_sat,
             'Q_GR2':Q_GR2_ref,
             'Q_GR3':Q_GR3_ref,
             'Q_GR4':Q_GR4_ref }])
    V_POI = pr.get_value('V_POI')
    
    # Controller
    error = V_POI_ref - V_POI
    xi = xi + T_s*error
    Q_POI_ref = K_p*error + K_i*xi
    Q_STATCOM_ref = Q_POI_ref
    if Q_STATCOM_ref >  Q_STATCOM_max: Q_STATCOM_ref =  Q_STATCOM_max
    if Q_STATCOM_ref < -Q_STATCOM_max: Q_STATCOM_ref = -Q_STATCOM_max
    Q_GR_ref = Q_POI_ref - Q_STATCOM_ref
    Q_GR1_ref = 0.3*Q_GR_ref
    Q_GR2_ref = 0.2*Q_GR_ref
    Q_GR3_ref = 0.3*Q_GR_ref
    Q_GR4_ref = 0.2*Q_GR_ref

    it += 1 

# Post-procesado
# ==============
pr.post()

fig1, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].step(times,out[:,1]/1e6,label='$Q_{POI}^\star$')
axes[0].step(times,out[:,2]/1e6,label='$Q_{STATCOM}^\star$')
axes[0].step(pr.T,pr.get_values('Q_GR1')/1e6,label='$Q_{GR1}^\star$')
axes[0].step(pr.T,pr.get_values('Q_GR2')/1e6,label='$Q_{GR2}^\star$')
axes[0].step(pr.T,pr.get_values('Q_GR3')/1e6,label='$Q_{GR3}^\star$')
axes[0].step(pr.T,pr.get_values('Q_GR4')/1e6,label='$Q_{GR4}^\star$')
axes[0].step(times,out[:,1]/1e6,label='$Q_{PMV}^\star$')

axes[1].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')

axes[0].set_ylabel('Reactive Power (Mvar)')
#axes[0].set_ylim([-5,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,10])

ax.set_xlabel('Time (s)')
fig1.align_ylabels()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### Control de tensión coordinado STATCOM+Generadores

# In[12]:


# Inicialización del sistema
# ==========================

Q_GR = 0e6
S_GR = 8e6

pr.initialize([{'P_GR1': 4e6, 'Q_GR1': Q_GR,
                'P_GR2': 6e6, 'Q_GR2': Q_GR,
                'P_GR3': 5e6, 'Q_GR3': Q_GR,
                'P_GR4':4.5e6, 'Q_GR4': Q_GR,
                'P_PMV': 0.0, 'Q_PMV': 0.0}], 1) 

K_g = 1/220e6 # gain of the plant

T_s = 0.1
times = np.arange(0,10,T_s)
N_t = len(times)
N_x = 2

# diseño del controlador
K_p = -1.0e7
K_i = 1/(T_s*K_g)

# condiciones iniciales
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
Q_GR1_ref = 0.0
Q_GR2_ref = 0.0
Q_GR3_ref = 0.0
Q_GR4_ref = 0.0
Q_STATCOM_ref = 0.0
V_POI_ref = pr.get_value('V_POI')

# vectores para guardar estados, etc.
out = np.zeros((N_t,2))
Q_GR_store = np.zeros((N_t,4))
Q_STATCOM_store = np.zeros((N_t,1))
it = 0

# Bucle de simulación
# ===================
for t in times:
    out[it,0] = xi
    out[it,1] = Q_POI_ref
    
    # perturbaciones o cambios de referencia
    if t>1.0:
        V_POI_ref = 1.05
        
    # planta
    for ig,q in zip(range(1,5),[Q_GR1_ref,Q_GR2_ref,Q_GR3_ref,Q_GR4_ref]):
        P_GR = pr.get_value(f'P_GR{ig}').real
        Q_GR_max = np.sqrt(S_GR**2 - P_GR**2)
        if q> Q_GR_max: q =  Q_GR_max
        if q<-Q_GR_max: q = -Q_GR_max
        pr.set_value(f'Q_GR{ig}',q)   
        Q_GR_store[it,ig-1] = q
        
    pr.run([{'t_end':t,'Q_PMV':Q_STATCOM_ref}])
    Q_STATCOM_store[it,0] = Q_STATCOM_ref
    V_POI = pr.get_value('V_POI')
    
    # controlador (PPC)
    error = V_POI_ref - V_POI
    xi = xi + T_s*error
    Q_POI_ref = K_p*error + K_i*xi
    p_poi = gt.get_line_s(pr,'POI','GRI').real
    q_gr_max = np.sqrt(S_GR**2 - (p_poi/4)**2)
    Q_GR = Q_POI_ref/4
    if Q_GR> q_gr_max: Q_GR =  q_gr_max
    if Q_GR<-q_gr_max: Q_GR = -q_gr_max
    Q_GR1_ref = Q_GR
    Q_GR2_ref = Q_GR
    Q_GR3_ref = Q_GR
    Q_GR4_ref = Q_GR
    
    Q_STATCOM_ref = Q_POI_ref - 4*Q_GR 
    
    it += 1 

# Post-procesado
# ==============
pr.post()
    
fig2, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].step(times,out[:,1]/1e6,label='$Q_{POI}^\star$')
for ig in range(4):
    axes[0].step(times,Q_GR_store[:,ig]/1e6,label=f'$Q_{{GR{ig+1}}}^\star$')
axes[0].step(times,Q_STATCOM_store[:,0]/1e6,label='$Q_{ST}^\star$')
axes[1].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')

axes[0].set_ylabel('Reactive Power (Mvar)')
axes[0].set_ylim([-5,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,10])

ax.set_xlabel('Time (s)')
fig2.align_ylabels()

