#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 16:35:34 2018

@author: ohm
"""

from bwhModel import bwhModel
import numpy as np



Ps_normal='auto/bwh_set2.hdf5'
Es_normal={'rhs':"oz",
        'n':(1024,),
        'l':(256.0,),
        'bc':"neumann",
        'it':"rk4",
        'dt':0.1,
        'analyze':False,
        'verbose':False,
        'setPDE':False}

def plotVars(Ps = Ps_normal,Es=Es_normal,
             Tstep=0.0001,Tstart=0,Tfinish=10,plot_p=False,
             initial_condition=[0.9,0.5,0.2]):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button#, RadioButtons
    m = bwhModel(Es=Es_normal,Ps=Ps,Vs=None)
    fig,ax=plt.subplots(3,1,sharex=True)
    ax[2].set_xlabel(r'$T \; [yr]$')
    ax[0].set_xlim([Tstart,Tfinish])
    ax[0].set_ylabel(r'$b$')
    ax[1].set_ylabel(r'$w$')
    ax[2].set_ylabel(r'$h$')
    plt.subplots_adjust(left=0.25, bottom=0.50)
    t,sol=m.ode_integrate(initial_state=initial_condition,step=Tstep*m.p['conv_T_to_t'],finish=Tfinish*m.p['conv_T_to_t'])
    idx = t.searchsorted(Tstart*m.p['conv_T_to_t'])
    Time=np.arange(Tstart,Tfinish+Tstep,Tstep)
    idx = np.clip(idx, 1, len(Time)-1)
    left = t[idx-1]
    right = t[idx]
    idx -= Tstart*m.p['conv_T_to_t'] - left < right - Tstart*m.p['conv_T_to_t']
    l1,=ax[0].plot(Time, sol[0][idx:],'g', lw=2)
    l2,=ax[1].plot(Time, sol[1][idx:],'b', lw=2)
    l3,=ax[2].plot(Time, sol[2][idx:],'c', lw=2)
    ax[1].axhline(y=m.p['w_wp'],linewidth=1, color = 'r',linestyle="--")
    if plot_p:
        pax, = ax.plot(t/m.p['conv_T_to_t'],m.p_t(t,m.p['p'],m.p['a'],m.p['omegaf']) , lw=1, color='b')
    axcolor = 'lightgoldenrodyellow'
    axchi  = plt.axes([0.25, 0.20, 0.65, 0.02])#, facecolor=axcolor)
    axtemp = plt.axes([0.25, 0.30, 0.65, 0.02])#, facecolor=axcolor)
    axp    = plt.axes([0.25, 0.25, 0.65, 0.02])#, facecolor=axcolor)
    axiota = plt.axes([0.25, 0.15, 0.65, 0.02])#, facecolor=axcolor)
    axomegaf = plt.axes([0.25, 0.10, 0.65, 0.02])#, facecolor=axcolor)
    schi  = Slider(axchi, r'$\chi$', 0., 1.0, valinit=m.p['chi'])
    stemp = Slider(axtemp, r'$T$', 1, 100.0, valinit=Tfinish)
    sp    = Slider(axp, r'$p$', 0.01, 10.0, valinit=m.p['p'])
    siota = Slider(axiota, r'$a$', 0.0, 1.0, valinit=m.p['a'])
    somegaf = Slider(axomegaf, r'$\omega_f$', 0.0, 50.0, valinit=m.p['omegaf'])

    def update(val):
        chi = schi.val
        Tfinish = stemp.val
        Time=np.arange(Tstart,Tfinish+Tstep,Tstep)
        a = siota.val
        ax[0].set_xlim([Tstart, Tfinish])
        p = sp.val
        omegaf = somegaf.val
        t,sol=m.ode_integrate(initial_state=initial_condition,p=p,chi=chi,a=a,step=Tstep*m.p['conv_T_to_t'],finish=Tfinish*m.p['conv_T_to_t'])
        idx = t.searchsorted(Tstart*m.p['conv_T_to_t'])
        Time=np.arange(Tstart,Tfinish+Tstep,Tstep)
        idx = np.clip(idx, 1, len(Time)-1)
        left = t[idx-1]
        right = t[idx]
        idx -= Tstart*m.p['conv_T_to_t'] - left < right - Tstart*m.p['conv_T_to_t']
        l1.set_ydata(sol[0])
        l1.set_xdata(Time)
        l2.set_ydata(sol[1])
        l2.set_xdata(Time)
        l3.set_ydata(sol[2])
        l3.set_xdata(Time)
        if plot_p:
            pax.set_ydata(m.p_t(t,p,a,omegaf))
            pax.set_xdata(t/m.p['conv_T_to_t'])
        fig.canvas.draw_idle()
    schi.on_changed(update)
    stemp.on_changed(update)
    sp.on_changed(update)
    siota.on_changed(update)
    somegaf.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button1 = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    def reset(event):
        schi.reset()
        stemp.reset()
        sp.reset()
        siota.reset()
    button1.on_clicked(reset)

    saveax = plt.axes([0.5, 0.025, 0.1, 0.04])
    button2 = Button(saveax, 'Save', color=axcolor, hovercolor='0.975')

    def save(event):
        chi = schi.val
        p = sp.val
        plt.savefig("results/bwh_set1_p_{:s}_chi_{:s}.png".format(str(int(p)),str(int(100.0*chi))))
    button2.on_clicked(save)

    plt.show()

if __name__ == '__main__':
    plotVars()