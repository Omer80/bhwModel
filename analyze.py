#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 07:08:40 2018

@author: ohm
"""
from bwhModel import bwhModel
import numpy as np

def simdraught(prec_i,prec_f,years,Ps,chi,
               n=(256,256),l=(128.0,128.0),
               Vs_initial="random",rhs="oz_EQK",
               bc="neumann",it="pseudo_spectral",
               first_time = 1000,tol=1.0e-8,
               fname="cont",verbose=True):
    import deepdish.io as dd
    Es={'rhs':rhs,'n':n,'l':l,'bc':bc,'it':it,
        'dt':0.1,'verbose':verbose,'analyze':False,'setPDE':True}
    if type(Vs_initial)==str:
        fname = fname+"_"+Vs_initial
    time_span = np.arange(0,years+1,1)
    m = bwhModel(Vs=Vs_initial,Es=Es,Ps=Ps)
    m.setup['verbose']=verbose
    # Converging on the first solution using integration and then root
    Vs_init = m.integrate(m.initial_state,check_convergence=True,
                          max_time=first_time,p=prec_i,chi=chi)
    Vs = Vs_init.copy()
    prec_gradient_down = np.linspace(prec_i,prec_f,len(time_span))
    b_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    w_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    h_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    for i,prec in enumerate(prec_gradient_down):
        print "Integration for p =",prec
        Vs_new=m.integrate(initial_state=Vs,max_time=m.p['conv_T_to_t'],
                           check_convergence=False,p=prec,chi=chi)
        b,w,h=m.split_state(Vs_new)
        b_sol[i]=b
        w_sol[i]=w
        h_sol[i]=h
        Vs = Vs_new
    dd.save(fname+".hdf5",{'p':prec_gradient_down,"T":time_span,
                           'b':b_sol,
                           'w':w_sol,
                           'h':h_sol},
            compression='blosc')

def followStable(prec_i,prec_min,prec_max,delta_prec,Ps,n=(128,128),l=(64.0,64.0),Vs_initial="random",root=False,track=False,max_time = 1000,tol=1.0e-8,fname="cont",verbose=True):
    import deepdish.io as dd
    Es={'rhs':"oz_EQK",'n':n,'l':l,'bc':"neumann",'it':"pseudo_spectral",
        'dt':0.1,'verbose':False,'analyze':False,'setPDE':True}
    if type(Vs_initial)==str:
        fname = fname+"_"+Vs_initial
    b_sol = []
    w_sol = []
    h_sol = []
    sol = []
    p_range = []
    conv_flag = []
    water_loss_ar = []
    m = bwhModel(Vs=Vs_initial,Es=Es,Ps=Ps)
#    m = afModel(Vs=Vs_initial,Es=Es,Ps=af_Ps)
    m.setup['verbose']=verbose
    # Converging on the first solution using integration and then root
    Vs_init = m.integrate_till_convergence(m.initial_state,prec_i)
    if track or root:
        Vs_init,converged = m.root(Vs_init,prec_i,fixiter=2000)
    Vs = Vs_init.copy()
    prec_gradient_up   = np.arange(prec_i,prec_max+delta_prec,delta_prec)
    prec_gradient_down = np.arange(prec_i-delta_prec,prec_min-delta_prec,-delta_prec)
    npoints = len(prec_gradient_up)+len(prec_gradient_down)
    # Setup netcdf4 continuation file
#    hn.setup_continuation(m,fname,npoints)
    for i,prec in enumerate(prec_gradient_up):
        if verbose:
            print "Process {0:.0f}%".format(float(i)/float(npoints) * 100)
            print "Finding solution for p=",prec
#        noutput = i+len(prec_gradient_down)
        if track:
            Vs_new,converged = m.root(Vs,prec,integrateiffail=True)
        elif root:
            Vs_new,converged = m.root(Vs,prec,fixiter=2000)
        else:
            Vs_new = m.integrate_till_convergence(m.initial_state,prec)
            converged = True
        b,w,h=m.splitfields(Vs_new)
        sol.append(Vs_new)
        b_sol.append(b)
        w_sol.append(w)
        h_sol.append(h)
        p_range.append(prec)
        conv_flag.append(converged)
        Vs = Vs_new
    Vs = Vs_init.copy() # Coming back to initial run
    for i,prec in enumerate(prec_gradient_down):
        if verbose:
            print "Process {0:.0f}%".format(float(i+len(prec_gradient_up))/float(npoints) * 100)
            print "Finding solution for p=",prec
        if track:
            Vs_new,converged = m.root(Vs,prec,integrateiffail=True)
        elif root:
            Vs_new,converged = m.root(Vs,prec,fixiter=2000)
        else:
            Vs_new = m.integrate_till_convergence(m.initial_state,prec)
            converged = True
        b,w,h=m.splitfields(Vs_new)
        sol.append(Vs_new)
        b_sol.append(b)
        w_sol.append(w)
        h_sol.append(h)
        p_range.append(prec)
        conv_flag.append(converged)
        Vs = Vs_new
    p_range=np.array(p_range)
    argssort = p_range.argsort()
    dd.save(fname+".hdf5",{'p':p_range[argssort],'n':n,'l':l,'conv_flag':np.array(conv_flag)[argssort],
            'b1_yield':np.array(b1_yield)[argssort],'b2_yield':np.array(b2_yield)[argssort],
            'Vs':np.array(sol)[argssort],'Ps':m.p,'Es':m.setup,
            'b1_sol':np.array(b1_sol)[argssort],'b2_sol':np.array(b2_sol)[argssort],
            'water_loss':np.array(water_loss_ar)[argssort]},
            compression='blosc')
