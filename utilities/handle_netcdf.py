# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:31:01 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
import time as t
import netCDF4
import numpy as np


def setup_simulation(fname,Ps,Es):
    """
    Opening an netCDF4 file with 3 groups:
    1 parameters - contains model's parameters set
    2 snapshots - will contain subgroups named after the continuation parameter
    value
    3 continuation - will contain datasets of array of stats values along the
    """
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
        for k,v in Ps.items():
            if k!='dimpar':
                setattr(parms,k,v)
        rootgrp.description = "Simulation dataset for tlmSpatial model."
        rootgrp.history = "Created " + t.ctime(t.time())
        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('time', None)
        time = rootgrp.createVariable('time', 'f8', ('time',),zlib=True)
        time.units = "year"
        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        x.units = "m"
        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        if len(Es['n']) == 1:
            print "Setting up 1D variables"
            for var in Es['Vs_symbols']:
                rootgrp.createVariable(str(var), 'f8', ('time', 'x',),zlib=True)
        elif len(Es['n']) == 2:
            print "Setting up 2D variables"
            rootgrp.createDimension("y", Es['n'][1])
            y = rootgrp.createVariable('y', 'f4', ('y',),zlib=True)
            y.units = "m"
            y[:] = np.linspace(0,Es['l'][1], Es['n'][1])
            for var in Es['Vs_symbols']:
                rootgrp.createVariable(str(var), 'f8', ('time', 'x', 'y',),zlib=True)
        print "Output: netCDF file was created: ", fname+".nc"
    
def save_sim_snapshot(fname,step,time,Vs,Es):
    """ Save snapshot of the four fields b1,b2,b3,s1,s2, together with the time
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['time'][step] = time
        if len(Es['n']) == 1:
            for i,var in enumerate(Es['Vs_symbols']):
                rootgrp[str(var)][step,:] = Vs[i]
        elif len(Es['n']) == 2:
            for i,var in enumerate(Es['Vs_symbols']):
                rootgrp[str(var)][step,:,:] = Vs[i]


def create_animation(fname,output=None,showtime=False):
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    # row and column sharing
    ims=[]
    with netCDF4.Dataset("%s.nc"%fname, 'r', format='NETCDF4') as rootgrp:
        nd = int(getattr(rootgrp['setup'],'nd'))
        t = rootgrp['time'][:]
        if nd == 1:
            fig, ax = plt.subplots(3, 1, sharex=True, sharey=True)
    #        axes = plt.gca()
            ax[0].set_ylim([-0.01,1.0])
            ax[1].set_ylim([-0.01,1.0])
            ax[2].set_ylim([-0.01,1.0])
            ax[0].set_ylabel(r'$b$', fontsize=25)
            ax[1].set_ylabel(r'$w$', fontsize=25)
            ax[2].set_ylabel(r'$h$', fontsize=25)
            ax[2].set_xlabel(r'$x$', fontsize=25)
            x  = rootgrp['x'][:]
            ax[0].set_xlim([x[0],x[-1]])
            b = rootgrp['b'][:,:]
            w = rootgrp['w'][:,:]
            h = rootgrp['h'][:,:]
            for i in xrange(len(t)):
                line1, = ax[0].plot(x,b[i],'g-')
                line2, = ax[1].plot(x,w[i],'g-')
                line3, = ax[2].plot(x,h[i],'g-')
                if showtime:
                    ax.set_title(r'$b$ at $t={:4.3f}$'.format(t[i]), fontsize=25)
                ims.append([line1,line2,line3])
        elif nd == 2:
            fig, ax = plt.subplots(1,2,sharex=True,sharey=True)
            fig.subplots_adjust(right=0.8)
            ax[0].set_aspect('equal', 'datalim')
            ax[0].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[0].set_title(r'$b$', fontsize=25)
            b = rootgrp['b'][:,:,:]
            ax[1].set_aspect('equal', 'datalim')
            ax[1].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[1].set_title(r'$w$', fontsize=25)
            w = rootgrp['w'][:,:,:]
            ax[2].set_aspect('equal', 'datalim')
            ax[2].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[2].set_title(r'$h$', fontsize=25)
            h = rootgrp['h'][:,:,:]
            for i in xrange(len(t)):
                im1, = ax.imshow(b[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                im2, = ax.imshow(w[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                im3, = ax.imshow(h[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                if showtime:
                    ax.set_title(r'$b$ at $t={:4.3f}$'.format(t[i]), fontsize=25)
                ims.append([im1,im2,im3])
                cbar_ax2 = fig.add_axes([0.85, 0.35, 0.05, 0.55])
                fig.colorbar(im1, cax=cbar_ax2)
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
    if output is None:
        ani.save('%s.mp4'%fname)
    elif type(output)==str:
        ani.save('%s.mp4'%output)
    print "Movie made!"

