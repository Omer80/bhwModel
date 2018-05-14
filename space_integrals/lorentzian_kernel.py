#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 17:13:10 2018

Giving me error when I try to use prange instead of Python's range

@author: ohm
"""
import numpy as np
from numba import njit, double, prange
import math

@njit(double[:](double[:],double[:],double,double,double,double,double), parallel=True)
def calc_gb_lorentz_1d(b,w,s0,e,lam,alpha,dx):
    n=b.shape[0]
    norm = 1.0/(np.pi*s0)
    gb = np.zeros((n))
    for i in prange(n):
        sigma=2.0*(s0*(1.0+e*b[i]))**2
        for ii in prange(n):
            gb[i]+=w[i]/((1.0+((((i-ii)*dx)**2)/sigma))**(1.0+alpha))
        gb[i]*=lam*norm
    return gb
@njit(double[:](double[:],double[:],double,double,double,double,double), parallel=True)
def calc_gw_lorentz_1d(b,w,s0,e,lam,alpha,dx):
    n=b.shape[0]
    norm = 1.0/(np.pi*s0)
    gw = np.zeros((n))
    for i in prange(n):
        for ii in prange(n):
            gw[i]+=b[ii]/((1.0+((((i-ii)*dx)**2)/(2.0*(s0*(1.0+e*b[ii]))**2)))**(1.0+alpha))
        gw[i]*=lam*norm
    return gw

#@njit(double[:, :](double[:, :],double[:,:],double,double,double,double), parallel=True)
#def calc_gb_gauss_2d(b,w,s0,e,lam,dx):
#    n,m=b.shape
#    norm = 1.0/(2*np.pi*s0**2)
#    gb = np.zeros((n,m))
#    for i in prange(n):
#        for j in prange(m):
#            sigma=2.0*(s0*(1.0+e*b[i,j]))**2
#            for ii in prange(n):
#                for jj in prange(m):
#                    gb[i,j]+=w[i,j]*math.exp(-(((i-ii)*dx)**2+((j-jj)*dx)**2)/sigma)
#            gb[i,j]*=lam*norm
#    return gb
#
#@njit(parallel=True)                 #needed, so it is parallelized
#def calc_gb_gauss_2d_ead(b,w,s0,e,lam,dx):
#    n,m=b.shape
#    norm = 1.0/(2*np.pi*s0**2)
#    gb = np.zeros((n,m))
#    vI=np.arange(n)
#    vJ=np.arange(m)
#    for i in prange(n):             #outer loop = explicit prange-loop
#        for j in range(m):
#            denom=2.0*(s0*(1.0+e*b[i,j]))**2
#            expII=np.zeros((n,))
#            expJJ=np.zeros((m,))
#            for k in range(n):
#                II=(i-vI[k])*dx
#                expII[k]=math.exp(-II*II/denom)
#
#            for k in range(m):
#                JJ=(j-vJ[k])*dx
#                expJJ[k]=w[i,j]*math.exp(-JJ*JJ/denom)
#            gb[i,j]=lam*norm*(expII.sum()*expJJ.sum())
#    return gb