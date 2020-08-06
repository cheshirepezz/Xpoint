#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:35:06 2019

author: Luca Pezzini
email: luca.pezzini@edu.unito.it
website: https://github.com/lucapezzini
license: BSD
"""

import matplotlib.pyplot as plt
import numpy as np

# Grid of x, y points
nx, ny = 20, 20
x = np.linspace(-1000, 1000, nx)
y = np.linspace(-1000, 1000, ny)
X, Y = np.meshgrid(x, y)
U = X/10
V = Y/10

plt.figure()
#plt.title('Arrows scale with plot width, not view')
M = np.hypot(U, V)


Q = plt.quiver(X, Y, U, V, M, units='width')
qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$1 B_{0}$', labelpos='E', coordinates='figure')

plt.scatter(X, Y, color = 'k', s = 5)
plt.axis(aspect = 'equal')
plt.xlabel('$x$', fontsize = 10)
plt.ylabel('$y$', fontsize = 10)
plt.savefig('Bfield.png', format='png', dpi=1200)
plt.show() 