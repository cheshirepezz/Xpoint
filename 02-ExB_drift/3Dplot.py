#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:28:40 2019

author: Luca Pezzini
email: luca.pezzini@edu.unito.it
website: https://github.com/lucapezzini
license: BSD
"""

import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

filename = "rk4.dat"
data = np.loadtxt(filename, delimiter = ' ')

x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

"""                                                                                                                                                    
Begin scaling                                                                                                                           
"""

x_scale = 1.5
y_scale = 1.3
z_scale = 1

scale = np.diag([x_scale, y_scale, z_scale, 1.0])
scale = scale*(1.0/scale.max())
scale[3,3] = 1.0

def short_proj():
  return np.dot(Axes3D.get_proj(ax), scale)

ax.get_proj=short_proj

"""                                                                                                                                                    
End scaling                                                                                                                           
"""
    
ax.scatter3D(x, y, z, label = 'RK4', c = 'blue', marker = '.')

#Shrink the number of coordinates in the graph
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
#zmin, zmax = ax.get_zlim()
ax.set_xticks(np.round(np.linspace(xmin, xmax, 5), 2))
ax.set_yticks(np.round(np.linspace(ymin, ymax, 5), 2))
#ax.set_zticks(np.round(np.linspace(zmin, zmax, 5), 2))

#Set the framework of the graph
ax.set_xlabel('$x$', fontsize = 10)
ax.set_ylabel('$y$', fontsize = 10)
ax.set_zlabel('$z$', fontsize = 10)
ax.set_title('ExB Drift', fontsize = 10)
ax.legend()

plt.savefig('RK4.png', format='png', dpi=1200)
plt.show() 