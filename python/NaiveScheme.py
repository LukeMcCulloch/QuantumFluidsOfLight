#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 13:10:53 2019

@author: lukemcculloch
"""

import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/17044052/mathplotlib-imshow-complex-2d-array
from colorsys import hls_to_rgb
#import matplotlib.colors.hsv_to_rgb as hls_to_rgp 

def colorize(z):
    n,m = z.shape
    c = np.zeros((n,m,3))
    c[np.isinf(z)] = (1.0, 1.0, 1.0)
    c[np.isnan(z)] = (0.5, 0.5, 0.5)

    idx = ~(np.isinf(z) + np.isnan(z))
    A = (np.angle(z[idx]) + np.pi) / (2*np.pi)
    A = (A + 0.5) % 1.0
    B = 1.0 - 1.0/(1.0+abs(z[idx])**0.3)
    c[idx] = [hls_to_rgb(a, b, 0.8) for a,b in zip(A,B)]
    return c

imag = np.complex(0.,1.)

class Grid(object):
    
    def __init__(self, hx=.01, hy=.01, Nx=100, Ny=100, lx=1., ly=1.):
        self.hx = hx
        self.hy = hy
        self.Nx = Nx
        self.Ny = Ny
        self.lx = lx
        self.ly = ly
        self.make_grid()
        
    def make_grid(self):
        """
        actual domain goes from 1-Nx, 1-Ny
        with boundary cells just beyond
        """
        self.grid = np.zeros((self.Nx+2,self.Ny+2),complex)
        self.gridupdate = np.zeros((self.Nx+2,self.Ny+2),complex)
        return
    

class TimeIntegrator(object):
    
    def __init__(self, grid, k):
        self.u = grid.grid
        self.u_update = grid.gridupdate
        self.Nx = grid.Nx
        self.Ny = grid.Ny
        self.dx = grid.hx
        self.dy = grid.hy
        self.k = k
        
    def shape(self):
        return np.shape(self.u)
        
    def DDx(self, i,j):
        u = self.u
        dx = self.dx
        return  (u[i-1,j] -2.*u[i,j] + u[i+1,j] )/dx
        
    def DDy(self, i,j):
        u = self.u
        dy = self.dy
        return  (u[i,j-1] -2.*u[i,j] + u[i,j+1] )/dy
    
    def IC(self, Q=1., i=None,j=None):
        if i is None:
            i = self.Nx/2 + 1
        if j is None:
            j = self.Ny/2 + 1
            
        self.u[i,j] = Q
        return 
    
    #def Dt(self):
    #    return 
    
    def solve(self, dt=.1, maxtime = 1.,
              ICfuncdict = None):
        #nstep = maxtime / dt
        if ICfuncdict is None:
            ICfuncdict = {}
            ICfuncdict[0] =  [self.IC, 1., 45, 45 ]
            ICfuncdict[1] =  [self.IC, 1., 55, 55 ]
            ICfuncdict[2] =  [self.IC, 1., 10, 55 ]
            ICfuncdict[3] =  [self.IC, 1., 55, 10 ]
            ICfuncdict[4] =  [self.IC, 1., 75, 75 ]
        
        t=0.
        sc = dt/(2.*imag*self.k)
        while( t<maxtime ):
            t += dt
            for i in range(1,self.Nx+1):
                for j in range(1,self.Ny+1):
                    self.u_update[i,j] = sc*( self.DDx(i,j) + self.DDy(i,j) ) + self.u[i,j]
                    for func in ICfuncdict:
                        clist = ICfuncdict[func]
                        clist[0]( Q = clist[1], i = clist[2], j = clist[3] )
            self.u[:,:] = self.u_update[:,:]
        return
    
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('colorMap')
        plt.imshow(colorize(self.u) )
        ax.set_aspect('equal')
        
        
        return
        
        
    
    
if __name__ == """__main__""":
    grid = Grid()
    self = grid
    
    ti = TimeIntegrator(grid, k=250.)
    ti.solve(maxtime=20.)
    ti.plot()