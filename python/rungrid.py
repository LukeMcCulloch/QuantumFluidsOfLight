# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:51:20 2019

@author: Luke.McCulloch
"""

# File rungrid.py: runs the multigrid programme.
import numpy as np
from grid import Grid
from smooth import get_lhs


n_grids=5
size_init=2**n_grids+1
size=size_init-1
h=1.0/size
foo=[] # Container list for grids


for k in range(n_grids):
    # Set up list of grids
    u=np.zeros((size+1,size+1),float)
    f=np.zeros_like(u)
    name='Level '+str(n_grids-1-k)
    temp=Grid(name,h,u,f)
    foo.append(temp)
    size/=2
    h*=2
    
    
for k in range(1,n_grids):
    # Set up coarser Grid links
    foo[k-1].co=foo[k]
    
    
    
# Check that the initial construction works
for k in range(n_grids):
    print foo[k]
    
# Set up data for the Numerical Recipes problem
u_init=np.zeros((size_init,size_init))
f_init=np.zeros_like(u_init)
f_init[size_init/2,size_init/2]=2.0
foo[0].u=u_init                             #trial solution
foo[0].f=f_init
foo[0].fmg_fas_v_cycle(1,1,1)

# As a check, get lhs for the final grid
lhs=get_lhs(foo[0].u,foo[0].h2)
print "max abs lhs = ", np.max(np.abs(lhs))


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib notebook



plt.ion()
fig1=plt.figure()
ax1=Axes3D(fig1)

xx,yy=np.mgrid[0:1:1j*size_init,0:1:1j*size_init]
ax1.plot_surface(xx,yy,1000*foo[0].u,rstride=1,cstride=1,alpha=0.2)
ax1.set_xlabel('x',style='italic')
ax1.set_ylabel('y',style='italic')
ax1.set_zlabel('1000*u',style='italic')

fig2=plt.figure()
ax2=Axes3D(fig2)
ax2.plot_surface(xx,yy,lhs,rstride=1,cstride=1,alpha=0.2)
ax1.set_xlabel('x',style='italic')
ax2.set_ylabel('y',style='italic')
ax2.set_zlabel('lhs',style='italic')

