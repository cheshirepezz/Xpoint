#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:28:40 2019

author: Luca Pezzini
email: luca.pezzini@edu.unito.it
website: https://github.com/lucapezzini
license: BSD
"""
#ffmpeg -framerate 30 -i RK4_%04d.png -c:v libx264 -r 30 out.mp4
#ffmpeg -framerate 24 -i RK4_%04d.png out.avi
import numpy as np
import matplotlib.pyplot as plt

#for i in range(1,10,1):
#filename = "rk4_1.dat"
filename = "boris_0030.dat"

data = np.loadtxt(filename, delimiter = ' ')
#data = np.loadtxt(filename%i, skiprows = 4, delimiter = ' ')

x = data[:, 1]
y = data[:, 2]
E = data[:, 7]

# Scatter the points, using size and color
plt.scatter(x, y, c = E, alpha = 0.5, cmap = 'jet', marker = '.')

#Shrink the number of coordinates in the graph
#xmin, xmax = plt.xlim()
#ymin, ymax = plt.ylim()
#Emin, Emax = plt.vlim()
xmin, xmax = - 1000, 1000
ymin, ymax = - 1000, 1000

plt.xticks(np.round(np.linspace(xmin, xmax, 5), 2))
plt.yticks(np.round(np.linspace(ymin, ymax, 5), 2))
#plt.vticks(np.round(np.linspace(vmin, vmax, 5), 2))


#Set the framework of the graph
plt.axis(aspect = 'equal')
plt.xlabel('$x$', fontsize = 10)
plt.ylabel('$y$', fontsize = 10)
plt.colorbar(label= 'E_k')

#plt.title('Time : %.4f'%(float(i)*0.15015), fontsize = 10)
   
#plt.savefig('RK4_1.eps', format='eps', dpi=1200)
#plt.savefig('RK4_0000.png', format='png', dpi=1200)
plt.show() 