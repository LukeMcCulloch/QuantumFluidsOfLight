#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 13:10:53 2019

@author: lukemcculloch
"""

import numpy as np


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
        self.grid = np.zeros((self.Nx+2,self.Ny+2),float)
        return
    

class TimeIntegrator(object):
    
    def __init__(self, grid):
        self.u = grid.grid
        self.Nx = grid.Nx
        self.Ny = grid.Ny
        self.dx = 
        
    def DDx(self, i,j):
        u = self.u
        return  (u[i-1,j] -2.*u[i,j] + u[i+1,j] )
    
    
if __name__ == """__main__""":
    grid = Grid()
    self = grid
    