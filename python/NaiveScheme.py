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
from matplotlib.colors import hsv_to_rgb

imag = np.complex(0.,1.)
pi = np.pi


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

def colorize2(z):
    r = np.abs(z)
    arg = np.angle(z) 

    h = (arg + pi)  / (2 * pi) + 0.5
    l = 1.0 - 1.0/(1.0 + r**0.3)
    s = 0.8

    c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.swapaxes(0,2) 
    return c

def Complex2HSV(z, rmin=-100000., rmax=100000., hue_start=90):
    # get amplidude of z and limit to [rmin, rmax]
    amp = np.abs(z)
    amp = np.where(amp < rmin, rmin, amp)
    amp = np.where(amp > rmax, rmax, amp)
    ph = np.angle(z, deg=1) + hue_start
    # HSV are values in range [0,1]
    h = (ph % 360) / 360
    s = 0.85 * np.ones_like(h)
    v = (amp -rmin) / (rmax - rmin)
    return hsv_to_rgb(np.dstack((h,s,v)))




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
    
    def IC(self, t, k, Q=1., i=None,j=None):
        omega = .05*pi
        if i is None:
            i = self.Nx/2 + 1
        if j is None:
            j = self.Ny/2 + 1
        
        #p = omega*imag*Q*np.exp(imag*k*t*2.*pi)
        p = omega*imag*Q
        #self.u[i,j] = p
        self.u_update[i-1,j] = p * np.exp(imag * (t * omega + pi*3.*.5)) #imag*pi
        self.u_update[i,j-1] = p * np.exp(imag * (t * omega + pi)) #imag*3.*pi*.5
        self.u_update[i+1,j] = p * np.exp(imag * (t * omega + pi*.5)) #
        self.u_update[i,j+1] = p * np.exp(imag * (t * omega )) #imag*pi*.5
        #self.u[i+1,j+1] = p
        #self.u[i-1,j-1] = p
        #self.u[i-1,j+1] = p
        #self.u[i+1,j-1] = p
        return 
    
    #def Dt(self):
    #    return 
    
    def solve(self, dt=.1, maxtime = 1.,
              ICfuncdict = None):
        #nstep = maxtime / dt
        if ICfuncdict is None:
            ICfuncdict = {}
            ICfuncdict[0] =  [self.IC, 10., self.Nx/2, self.Ny/2-self.Ny/10 ]
            ICfuncdict[1] =  [self.IC, -10., self.Nx/2, self.Ny/2+self.Ny/10 ]
            
#            ICfuncdict[2] =  [self.IC, 10., self.Nx/4, self.Ny/2-self.Ny/10 ]
#            ICfuncdict[3] =  [self.IC, -10., self.Nx/4, self.Ny/2+self.Ny/10 ]
#            ICfuncdict[4] =  [self.IC, 10., 3*self.Nx/4, self.Ny/2-self.Ny/10 ]
#            ICfuncdict[5] =  [self.IC, -10., 3*self.Nx/4, self.Ny/2+self.Ny/10 ]
        
        t=0.
        for func in ICfuncdict:
            clist = ICfuncdict[func]
            clist[0]( t, k=self.k, Q = clist[1], i = clist[2], j = clist[3] )
            
        sc = dt/(2.*imag*self.k)
        while( t<maxtime ):
            t += dt
            for i in range(1,self.Nx+1):
                for j in range(1,self.Ny+1):
                    self.u_update[i,j] = sc*( self.DDx(i,j) + self.DDy(i,j) ) + self.u[i,j]
            for func in ICfuncdict:
                clist = ICfuncdict[func]
                clist[0]( t, k=self.k, Q = clist[1], i = clist[2], j = clist[3] )
            self.u[:,:] = self.u_update[:,:]
        return
    
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Amplitude and Phase')
        plt.imshow(colorize(self.u) , interpolation='quadric')
        ax.set_aspect('equal')
        return
        
    def plot_amp(self):
        z = self.u
        r = np.abs(z)
        #arg = np.angle(z) 
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('amplitude')
        plt.imshow(r , interpolation='none',extent=(-5,5,-5,5))
        ax.set_aspect('equal')
        return
        
    
    def plot_phase(self):
        z = self.u
        #r = np.abs(z)
        arg = np.angle(z) 
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('phase')
        plt.imshow(arg , interpolation='none',extent=(-5,5,-5,5))
        ax.set_aspect('equal')
        return
        
        
    
    
if __name__ == """__main__""":
    grid = Grid(Nx=150, Ny=150)
    self = grid
    
    ti = TimeIntegrator(grid, k=250.)
    ti.solve(dt=.1, maxtime=20.)
    ti.plot()
    #ti.plot_amp()
    #ti.plot_phase()
    
    
#N = 100
#A = np.zeros((N,N),dtype='complex')
#axis_x = np.linspace(-5,5,N)
#axis_y = np.linspace(-5,5,N)
#X,Y = np.meshgrid(axis_x,axis_y)
#Z = X + Y*1j
#
#A = 1/(Z+1j)**2 + 1/(Z-2)**2
#
## Plot the array "A" using colorize
#import pylab as plt
#plt.imshow(colorize(A), interpolation='none',extent=(-5,5,-5,5))
#plt.show()