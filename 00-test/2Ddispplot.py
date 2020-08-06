#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:28:40 2019

author: Luca Pezzini
email: luca.pezzini@edu.unito.it
website: https://github.com/lucapezzini
license: BSD
"""
import matplotlib.pyplot as plt

filename1 = "rk4_err.dat"
filename2 = "boris_err.dat"

data1 = np.loadtxt(filename1, delimiter = ' ')
data2 = np.loadtxt(filename2, delimiter = ' ')

x1 = data1[:, 0]
y1 = data1[:, 2]
x2 = data2[:, 0]
y2 = data2[:, 2]

# Scatter the points, using size and color
plt.loglog(x1, y1, c = 'red', marker = '.', label='RK4')
plt.loglog(x2, y2, c = 'blue', marker = '.', label='Boris')
#Set the framework of the graph
plt.axis(aspect = 'equal')
plt.xlabel('$log(h)$', fontsize = 10)
plt.ylabel('$log(err)$', fontsize = 10)
plt.legend()
#plt.savefig('relx.png', format='png', dpi=1200)
#plt.savefig('relv.png', format='png', dpi=1200)
plt.show() 