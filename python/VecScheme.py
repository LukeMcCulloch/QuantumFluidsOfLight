#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 13:10:53 2019

@author: lukemcculloch
"""

import numpy as np
from numpy import exp
#import scipy as sp
import scipy.special as special
Hermite = special.hermite
#laplace = sp.ndimage.filters.laplace
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
    
    #A = (np.angle(z[idx]) ) 
    B = 1.0 - 1.0/(1.0+abs(z[idx])**0.1)
    #B = abs(z[idx])
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



class Forcing(object):
    
    def __init__(self, location, N, size):
        self.N = N
        self.size = size
        self.location = location
        
    
    
    def gaussian(Lx=10, Ly=10,c=1.0):
        """
        Initial Gaussian bell in the middle of the domain.
        """
    
    
        def I(x, y):
            """Gaussian peak at (Lx/2, Ly/2)."""
            return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
    
        return I
    
    
    
    def gaussian2(self, Lx=10, Ly=10,c=1.0):
        """
        Initial Gaussian bell in the domain.
        """
        
    
    
        def I(x, y, w):
            """Gaussian peak at (Lx/2, Ly/2)."""
            return exp( -(x**2+y**2) /  w**2 )
        return I
    
    def GaussHermite( self, n, m, A, w0, Field ):
        N = self.N
        size = self.size
    
        sqrt2w0 = np.sqrt(2.0)/w0
        w02     = w0*w0
        n2      = N/2
        dx      = size/N
    
        for i in range(0,N):
            x = (i-n2)*dx
            x2 = x*x
            sqrt2xw0 = sqrt2w0*x
            
            for j in range(0,N):
                y = (j-n2)*dx
                y2 = y*y
                sqrt2yw0 = sqrt2w0*y
                Field[i,j ] = (A*exp(-(x2+y2)/w02)*Hermite(m,sqrt2xw0)*Hermite(n,sqrt2yw0) , 0.0)
        
        
        return Field
    
    





class Grid(object):
    
    def __init__(self, 
                 hx=.01, hy=.01, 
                 Nx=100, Ny=100, 
                 lx=1.,  ly=1.):
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
        self.shape = np.shape(self.u)
        self.u_update = grid.gridupdate
        self.Nx = grid.Nx
        self.Ny = grid.Ny
        self.dx = grid.hx
        self.dy = grid.hy
        self.k = k
        
    def shape(self):
        return np.shape(self.u)
    
        
    
    def DDx(self):
        nx = self.shape[0]
        ny = self.shape[1]
        u = self.u
        return   u[:-2,1:-1] -2.*u[1:-1,1:-1] + u[2:,1:-1]  


        
    def DDy(self):
        nx = self.shape[0]
        ny = self.shape[1]
        u = self.u
        return  u[1:-1,:-2] -2.*u[1:-1,1:-1] + u[1:-1,2:] 

    
    
    def IC(self, t, k, Q=1., omega=.05*pi, i=None,j=None):
        #omega = 0.15*pi
        if i is None:
            i = self.Nx/2 + 1
        if j is None:
            j = self.Ny/2 + 1
        
        #p = omega*imag*Q*np.exp(imag*k*t*2.*pi)
        p = Q #
        #p = omega*imag*Q
        #self.u_update[i,j] = p
        self.u_update[i-1,j] = p * np.exp(imag * (t * omega + pi*3.*.5)) #imag*pi
        self.u_update[i,j-1] = p * np.exp(imag * (t * omega + pi)) #imag*3.*pi*.5
        self.u_update[i+1,j] = p * np.exp(imag * (t * omega + pi*.5)) #
        self.u_update[i,j+1] = p * np.exp(imag * (t * omega )) #imag*pi*.5
        
#        self.u[i+1,j+1] = p * np.exp(imag * (t * omega + .25*pi))
#        self.u[i-1,j-1] = p * np.exp(imag * (t * omega + .75*pi))
#        self.u[i-1,j+1] = p * np.exp(imag * (t * omega + 1.25*pi))
#        self.u[i+1,j-1] = p * np.exp(imag * (t * omega + 1.75*pi))
        return 
    
    def iIC(self, t, k, Q=1., omega = .05*pi, i=None,j=None):
        #omega = .1*pi
        if i is None:
            i = self.Nx/2 + 1
        if j is None:
            j = self.Ny/2 + 1
        
        #p = omega*imag*Q*np.exp(imag*k*t*2.*pi)
        p = Q #
        #p = omega*imag*Q
        #self.u[i,j] = p
        self.u_update[i-1,j] = p * np.exp(imag * (t * omega - pi*3.*.5)) #imag*pi
        self.u_update[i,j-1] = p * np.exp(imag * (t * omega - pi)) #imag*3.*pi*.5
        self.u_update[i+1,j] = p * np.exp(imag * (t * omega - pi*.5)) #
        self.u_update[i,j+1] = p * np.exp(imag * (t * omega )) #imag*pi*.5
        
        #        self.u[i+1,j+1] = p * np.exp(imag * (-t * omega - .25*pi))
        #        self.u[i-1,j-1] = p * np.exp(imag * (-t * omega - .75*pi))
        #        self.u[i-1,j+1] = p * np.exp(imag * (-t * omega - 1.25*pi))
        #        self.u[i+1,j-1] = p * np.exp(imag * (-t * omega - 1.75*pi))
        return 
    
    #def Dt(self):
    #    return 
    
    def solve(self, dt=.1, maxtime = 1.,
              ICfuncdict = None, w=1.9):
        
        t=0.
        for func in ICfuncdict:
            clist = ICfuncdict[func]
            clist[0]( t, k=self.k, Q = clist[1], i = clist[2], j = clist[3] )
        #self.u[:,:] = np.copy(self.u_update[:,:])
            
        sc = dt/(2.*imag*self.k)
        
        invdy = Cy2 = 1./self.dy
        invdx = Cx2 = 1./self.dx
        
        while( t<maxtime ):
            u_n = np.copy(self.u)
            
            t += dt
            
               
            self.u_update[1:-1,1:-1] = u_n[1:-1,1:-1] + \
                                             sc* (  Cx2*self.DDx() + Cy2*self.DDy() ) 
            self.u_update[1:-1,1:-1] =  w*self.u_update[1:-1,1:-1] + \
                                                (1.-w)*self.u[1:-1,1:-1]
            

            for func in ICfuncdict:
                clist = ICfuncdict[func]
                clist[0]( t, k=self.k, Q = clist[1], i = clist[2], j = clist[3] )
            
                
            self.u[:,:] = self.u_update[:,:]
        return
    
    def plot(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Amplitude and Phase')
        #plt.imshow(colorize(self.u) , interpolation='none')
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
    
    grid = Grid(Nx=300, Ny=300)
    self = grid
    
    ti = TimeIntegrator(grid, k=1000.)
    
    ICfuncdict = {}
    #ICfuncdict[0] =  [ti.IC,  10., ti.Nx/2, ti.Ny/2]#
    #ICfuncdict[0] =  [ti.IC,  10., ti.Nx/2, ti.Ny/2-ti.Ny/10 ]
    ICfuncdict[0] =  [ti.IC,  10., ti.Nx/2, ti.Ny/2-ti.Ny/10 ]
    ICfuncdict[1] =  [ti.iIC, -10., ti.Nx/2, ti.Ny/2+ti.Ny/10 ]
    
    ICfuncdict[2] =  [ti.IC, 10., ti.Nx/3, ti.Ny/2]
    #ICfuncdict[3] =  [ti.IC, -10., ti.Nx/4, ti.Ny/2+ti.Ny/10 ]
    #ICfuncdict[4] =  [ti.IC, 10., 3*ti.Nx/4, ti.Ny/2-ti.Ny/10 ]
    #ICfuncdict[5] =  [ti.IC, -10., 3*ti.Nx/4, ti.Ny/2+ti.Ny/10 ]
    
    ti.solve(dt=.1, 
             maxtime=200.,
             ICfuncdict=ICfuncdict,
             w=1.99)
    ti.plot()
    #ti.plot_amp()
    #ti.plot_phase()
    self = ti
    
"""
N = 100
A = np.zeros((N,N),dtype='complex')
axis_x = np.linspace(-5,5,N)
axis_y = np.linspace(-5,5,N)
X,Y = np.meshgrid(axis_x,axis_y)
Z = X + Y*1j

A = 1/(Z+1j)**2 + 1/(Z-2)**2

# Plot the array "A" using colorize
import pylab as plt
plt.imshow(colorize(A), interpolation='quadric',extent=(-5,5,-5,5))
plt.show()
#"""