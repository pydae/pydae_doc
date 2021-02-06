#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from pydgrid.plot_bokeh import plot_results
import sympy as sym
import pydae.build as db
from pydae.grid_urisi import unb_ri_si
import json


# In[2]:


data = {
        "buses":[
                 {"bus": "B1",  "pos_x":   0, "pos_y":  0, "units": "m", "U_kV":0.4},
                 {"bus": "B2",  "pos_x":  20, "pos_y":  0, "units": "m", "U_kV":0.4},
                 {"bus": "B3",  "pos_x": 120, "pos_y":  0, "units": "m", "U_kV":0.4},
                 {"bus": "B4",  "pos_x": 140, "pos_y":  0, "units": "m", "U_kV":0.4}
                ],
        "grid_formers":[
                        {"bus": "B1",
                        "bus_nodes": [1, 2, 3], "deg": [0, -120, -240],
                        "kV": [0.231, 0.231, 0.231]},
                        {"bus": "B4",
                        "bus_nodes": [1, 2, 3], "deg": [0, -120, -240],
                        "kV": [0.231, 0.231, 0.231]}
                       ],
        "lines":[
                 {"bus_j": "B1",  "bus_k": "B2",  "code": "lv_cu_150", "m":  20.0},
                 {"bus_j": "B2",  "bus_k": "B3",  "code": "lv_cu_150", "m": 100.0},
                 {"bus_j": "B3",  "bus_k": "B4",  "code": "lv_cu_150", "m":  20.0},
                ],
        "loads":[
                 {"bus": "B2" , "kVA": [30.0,30.0,30.0], "pf":[ 1]*3,"type":"3P+N"},
                 {"bus": "B3" , "kVA": [10.0,10.0,70.0], "pf":[ 1]*3,"type":"3P+N"}
                ],
        "shunts":[
                 {"bus": "B1" , "R": 0.001, "X": 0.0, "bus_nodes": [4,0]},
                 {"bus": "B4" , "R": 0.001, "X": 0.0, "bus_nodes": [4,0]}
                 ],
        "line_codes":
            {"lv_cu_150":  {"Rph":0.167,"Xph":0.08, "Rn":0.167, "Xn": 0.08}
            }
       }


# In[3]:


grid_dae = unb_ri_si(data)


# In[4]:


params_dict  = grid_dae['params']
f_list = grid_dae['f']
x_list = grid_dae['x']
g_list = grid_dae['g']
y_list = grid_dae['y']
u_dict = grid_dae['u']
a = sym.Symbol('a')
h_dict = grid_dae['h_v_m_dict']

sys_dict = {'name':'grid_4bus4wire',
           'params_dict':params_dict,
           'f_list':f_list,
           'g_list':g_list,
           'x_list':x_list,
           'y_ini_list':y_list,
           'y_run_list':y_list,
           'u_run_dict':u_dict,
           'u_ini_dict':u_dict,
           'h_dict':h_dict
           }

db.system(sys_dict)
db.sys2num(sys_dict)

data = json.dumps(grid_dict['xy_0_dict'], indent=4)
fobj = open("xy_0_dict.json","w")
fobj.write(data)
fobj.close()


# In[ ]:


Y_ii = grid_1.Y_ii.toarray()
Y_vv = grid_1.Y_vv
Y_vi = grid_1.Y_vi
inv_Y_ii = np.linalg.inv(Y_ii)
N_nz_nodes = grid_1.params_pf[0].N_nz_nodes
N_v = grid_1.params_pf[0].N_nodes_v
nodes_list = grid_1.nodes
Y_primitive = grid_1.Y_primitive_sp.toarray() 
A_conect = grid_1.A_sp.toarray()
node_sorter  = grid_1.node_sorter
N_v = grid_1.N_nodes_v

np.savez('matrices',Y_primitive=Y_primitive,A_conect=A_conect,nodes_list=nodes_list,
         node_sorter=node_sorter,N_v=N_v, Y_vv=Y_vv, Y_vi=Y_vi)


with open("grid_data.json", "w") as fobj:
    json.dump(grid_1.data, fobj, indent=4, sort_keys=True)


# In[ ]:





# In[ ]:




