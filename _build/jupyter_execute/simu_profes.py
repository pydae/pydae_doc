#!/usr/bin/env python
# coding: utf-8

# # Ejemplo de simulación

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from pydae.systems.ier.proyecto import proyecto_class
import pydae.grid_tools as gt
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'svg'")


# Primero se instancia el simulador:

# In[2]:


pr = proyecto_class()


# Los parametros de las líneas y transformadores deben ser adaptados, para ello se puede utilizar la herramienta `change_line` del modulo `pydae.grid_tools`:

# In[3]:


gt.change_line(pr,'GR4','GR3',R_km=8.0,X_km=4.0,km=0.5)
gt.change_line(pr,'GR2','GR1',R_km=8.0,X_km=4.0,km=0.5)
gt.change_line(pr,'GR3','PMV',R_km=4.0,X_km=2.0,km=0.5)
gt.change_line(pr,'GR1','PMV',R_km=4.0,X_km=2.0,km=0.5)
gt.change_line(pr,'PMV','POI',R_pu=0.01,X_pu=0.05,S_mva=52.0)
gt.change_line(pr,'POI','GRI',R_km=0.168495,X_km=0.04,km=26)


# Una vez que el sistema tiene los parametros deseados se lo puede inicializar.

# In[4]:


P_GR = 5e6
Q_GR = 5e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR}], 1) 


# Para comprobar que la inicialización fue exitosa se pueden consultar los valores de las variables algebráicas:

# In[5]:


pr.report_y()


# In[6]:


P_GR = 5e6
Q_GR = 0e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR}], 1) 
pr.initialize([{'Q_PMV':0.0}],1)
print(f"V_POI antes del cambio: {pr.get_value('V_POI')}")


# Estimación de la potencia de cortocicuito
# 

# In[7]:


Q_PMV_1 = 0.0
Q_PMV_2 = 10.0e6
pr.initialize([{'Q_PMV':Q_PMV_1}],1)
V_POI_1 = pr.get_value('V_POI')
pr.initialize([{'Q_PMV':Q_PMV_2}],1)
V_POI_2 = pr.get_value('V_POI')

S_cc = (Q_PMV_2 - Q_PMV_1)/(V_POI_2 - V_POI_1)
S_cc/1e6


# ## Control de tensión del POI
# 
# ### Control de tensión utilizando sólo el STATCOM

# In[8]:


P_GR = 5e6
Q_GR = 0e6

pr.initialize([{'P_GR1': P_GR, 'Q_GR1': Q_GR,
                'P_GR2': P_GR, 'Q_GR2': Q_GR,
                'P_GR3': P_GR, 'Q_GR3': Q_GR,
                'P_GR4': P_GR, 'Q_GR4': Q_GR,
                'P_PMV': 0.0, 'Q_PMV': 0.0}], 1) 

K_g = 1/300e6 # gain of the plant

T_s = 0.1
times = np.arange(0,10,T_s)
N_t = len(times)
N_x = 2

# controller
K_p = 0.0
K_i = 1/(T_s*K_g)

# initial conditions
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
V_POI_ref = pr.get_value('V_POI')

X = np.zeros((N_t,N_x))
it = 0
for t in times:
    X[it,0] = xi
    X[it,1] = Q_POI_ref
    
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

pr.post()

fig1, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].step(times,X[:,1]/1e6,label='$Q_{POI}^\star$')
axes[0].step(times,X[:,1]/1e6,label='$Q_{PMV}^\star$')

axes[1].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')

axes[0].set_ylabel('Reactive Power (Mvar)')
axes[0].set_ylim([0,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,10])

ax.set_xlabel('Time (s)')
fig1.align_ylabels()


# ### Control de tensión utilizando coordinando generadores y STATCOM

# In[9]:


P_GR = 7e6
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

# controller
K_p = -1.0e7
K_i = 1/(T_s*K_g)

# initial conditions
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
Q_GR1_ref = 0.0
Q_GR2_ref = 0.0
Q_GR3_ref = 0.0
Q_GR4_ref = 0.0
Q_STATCOM_ref = 0.0
V_POI_ref = pr.get_value('V_POI')
Q_GR_store = np.zeros((N_t,4))
Q_STATCOM_store = np.zeros((N_t,1))
X = np.zeros((N_t,N_x))
it = 0
for t in times:
    X[it,0] = xi
    X[it,1] = Q_POI_ref
    
    if t>1.0:
        V_POI_ref = 1.05
        
    # Plant
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
    
    # Controller
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

pr.post()
    
fig2, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].step(times,X[:,1]/1e6,label='$Q_{POI}^\star$')
for ig in range(4):
    axes[0].step(times,Q_GR_store[:,ig]/1e6,label=f'$Q_{{GR{ig+1}}}^\star$')
axes[0].step(times,Q_STATCOM_store[:,0]/1e6,label='$Q_{ST}^\star$')
axes[1].plot(pr.T,pr.get_values('V_POI'),label='$V_{POI}$')

axes[0].set_ylabel('Reactive Power (Mvar)')
axes[0].set_ylim([0,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,10])

ax.set_xlabel('Time (s)')
fig2.align_ylabels()


# In[10]:


4+6+5+4.5


# In[11]:


P_GR = 7e6
Q_GR = 0e6
S_GR = 8e6
P_GR_n = 6e6
P_GRI = -100e6
pr.initialize([{'P_GR1': 4e6, 'Q_GR1': Q_GR,
                'P_GR2': 6e6, 'Q_GR2': Q_GR,
                'P_GR3': 5e6, 'Q_GR3': Q_GR,
                'P_GR4':4.5e6, 'Q_GR4': Q_GR,
                'P_PMV': 0.0, 'Q_PMV': 0.0,
                'P_GRI': P_GRI, 'Q_GRI': 50e6,'S_n_GRI':200.0e6}], 1) 

K_g = 1/220e6 # gain of the plant

T_s = 0.1
times = np.arange(0,30,T_s)
N_t = len(times)
N_x = 2

# controller
K_p = -1.0e7
K_i = 1/(T_s*K_g)

# initial conditions
x = 0.0
xi = 0.0
Q_POI_ref = 0.0
Q_GR1_ref = 0.0
Q_GR2_ref = 0.0
Q_GR3_ref = 0.0
Q_GR4_ref = 0.0
Q_STATCOM_ref = 0.0
P_BESS_ref = 0.0
V_POI_ref = pr.get_value('V_POI')
Q_GR_store = np.zeros((N_t,4))
Q_STATCOM_store = np.zeros((N_t,1))
X = np.zeros((N_t,N_x))
it = 0
for t in times:
    X[it,0] = xi
    X[it,1] = Q_POI_ref
    
    if t>100.0:
        V_POI_ref = 1.05
    if t>1.0:
        P_GRI = -150e6 
        
    # Plant
    for ig,q in zip(range(1,5),[Q_GR1_ref,Q_GR2_ref,Q_GR3_ref,Q_GR4_ref]):
        P_GR = pr.get_value(f'P_GR{ig}').real
        Q_GR_max = np.sqrt(S_GR**2 - P_GR**2)
        if q> Q_GR_max: q =  Q_GR_max
        if q<-Q_GR_max: q = -Q_GR_max
        pr.set_value(f'Q_GR{ig}',q)   
        Q_GR_store[it,ig-1] = q
    pr.set_value(f'P_GRI',P_GRI)    
    Q_STATCOM_store[it,0] = Q_STATCOM_ref
    V_POI = pr.get_value('V_POI')
    
    
    # run
    pr.run([{'t_end':t,'Q_PMV':Q_STATCOM_ref, 'P_PMV':P_BESS_ref}])
    
    
    # Controller
    error = V_POI_ref - V_POI
    xi = xi + T_s*error
    Q_POI_ref = K_p*error + K_i*xi
    p_poi = gt.get_line_s(pr,'POI','GRI').real
    if p_poi/4 < S_GR: 
        q_gr_max = np.sqrt(S_GR**2 - (p_poi/4)**2)
    else: q_gr_max = np.sqrt(S_GR**2 - P_GR_n**2) 
    Q_GR = Q_POI_ref/4
    if Q_GR> q_gr_max: Q_GR =  q_gr_max
    if Q_GR<-q_gr_max: Q_GR = -q_gr_max
    Q_GR1_ref = Q_GR
    Q_GR2_ref = Q_GR
    Q_GR3_ref = Q_GR
    Q_GR4_ref = Q_GR
    Q_STATCOM_ref = Q_POI_ref - 4*Q_GR 
    
    omega = pr.get_value('omega_GRI')
    p_f = 1/0.05*(1 - omega)
    P_BESS_ref = p_f*50e6
    it += 1 

pr.post()
    
fig2, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 6),sharex=True)

axes[0].plot(pr.T,pr.get_values('omega_GRI')*50,label='freq. (Hz)')

for ig in range(4):
    axes[1].plot(pr.T,pr.get_values(f'P_GR{ig+1}')/1e6,label=f'$P_{{GR{ig+1}}}$')
axes[1].plot(pr.T,pr.get_values(f'P_PMV')/1e6,label=f'$P_{{PMV}}$')
#axes[0].set_ylabel('Reactive Power (Mvar)')
#axes[0].set_ylim([-5,35])

axes[1].set_ylabel('Voltage (p.u.)')
for ax in axes:
    ax.grid()
    ax.legend(loc='best',ncol=3)
    ax.set_xlim([0,30])

ax.set_xlabel('Time (s)')
fig2.align_ylabels()


# In[48]:


pr.report_x()


# In[60]:


pr.get_value('v_f_GRI')


# In[45]:


S_GR**2 - (p_poi/4)**2


# In[46]:


p_poi/1e6


# In[37]:


43750000000000.0/1e6


# In[38]:


S_GR


# In[39]:


P_GR


# In[41]:


np.sqrt(S_GR**2 - P_GR**2)


# In[ ]:





# In[21]:


P_GR = 5e6
P_nom = 4*P_GR # potencia nominal de nuestro parque
print(f'P_nom = {P_nom/1e6} MW')


# In[12]:


E_alcolea  = 3_300_000e3 # Wh
P_nom_alcolea =  2_160e3 # W


# In[23]:


E_planta = P_nom/P_nom_alcolea*E_alcolea
E_planta_MWh = E_planta/1e6
print(f'E_planta_MWh = {E_planta_MWh:0.0f} MWh')
E_planta_MWh


# In[25]:


euro_MWh = 53.41
euro_año = E_planta_MWh*euro_MWh
print(f'euro_año = {euro_año/1e6:0.1f} M€')


# In[31]:


euros_sin_bateria = 0.1*euro_año 

print(f'euros sin bateria = {euros_sin_bateria/1e3:0.2f} k€')
print(f'total euros sin bateria = {20*euros_sin_bateria/1e6:0.2f} M€')


# In[40]:


p_bateria = 0.1*P_nom
e_bateria_J = 2*p_bateria * (10*60)
e_bateria_Wh = e_bateria_J/3600
print(f'energía batería = {e_bateria_Wh/1e3:0.2f} kWh')

euro_Wh_batería = 1.0
costo_batería = e_bateria_Wh*euro_Wh_batería
print(f'costo batería = {costo_batería/1e3:0.2f} k€')
print(f'costo total con batería sin curtailment = {2*costo_batería/1e6:0.2f} M€')


# In[42]:


p_bateria = 0.1*P_nom
e_bateria_J = p_bateria * (10*60)
e_bateria_Wh = e_bateria_J/3600
print(f'energía batería = {e_bateria_Wh/1e3:0.2f} kWh')

euro_Wh_batería = 1.0
costo_batería = e_bateria_Wh*euro_Wh_batería
print(f'costo batería = {costo_batería/1e3:0.2f} k€')
print(f'costo total con batería con curtailment = {2*costo_batería/1e6:0.2f} M€')


# In[20]:


costo_curt = 0.1*euro_año
costo_curt


# In[ ]:




