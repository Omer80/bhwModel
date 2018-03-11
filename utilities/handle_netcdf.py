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
        rootgrp.createVariable('b1_coverage', 'f8', ('time',),zlib=True)
        rootgrp.createVariable('b2_coverage', 'f8', ('time',),zlib=True)
        rootgrp.createVariable('b3_coverage', 'f8', ('time',),zlib=True)
        if len(Es['n']) == 1:
            print "Setting up 1D variables"
            rootgrp.createVariable('b1', 'f8', ('time', 'x',),zlib=True)
            rootgrp.createVariable('b2', 'f8', ('time', 'x', ),zlib=True)
            rootgrp.createVariable('b3', 'f8', ('time', 'x', ),zlib=True)
            rootgrp.createVariable('s1', 'f8', ('time', 'x', ),zlib=True)
            rootgrp.createVariable('s2', 'f8', ('time', 'x', ),zlib=True)
        elif len(Es['n']) == 2:
            print "Setting up 2D variables"
            rootgrp.createDimension("y", Es['n'][1])
            y = rootgrp.createVariable('y', 'f4', ('y',),zlib=True)
            y.units = "m"
            y[:] = np.linspace(0,Es['l'][1], Es['n'][1])
            rootgrp.createVariable('b1', 'f8', ('time', 'x', 'y',),zlib=True)
            rootgrp.createVariable('b2', 'f8', ('time', 'x', 'y',),zlib=True)
            rootgrp.createVariable('b3', 'f8', ('time', 'x', 'y',),zlib=True)
            rootgrp.createVariable('s1', 'f8', ('time', 'x', 'y',),zlib=True)
            rootgrp.createVariable('s2', 'f8', ('time', 'x', 'y',),zlib=True)
        print "Output: netCDF file was created: ", fname+".nc"
def setup_p_a_scan(fname,Ps,Es,t_range,p_range,a_range):
    """
    Opening an netCDF4 file with 3 groups:
    1 parameters - contains model's parameters set
    2 snapshots - will contain subgroups named after the continuation parameter
    value
    3 continuation - will contain datasets of array of stats values along the
    """
    assert (len(Es['n']) == 1),"Able to save p to a scan only for 1D simulations"
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
        rootgrp.description = "p to theta grid scan dataset for onfc2s model."
        rootgrp.history = "Created " + t.ctime(t.time())
        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('time', len(t_range))
        rootgrp.createDimension('p', len(p_range))
        rootgrp.createDimension('a', len(a_range))
        time = rootgrp.createVariable('time', 'f4', ('time',),zlib=True)
        time.units = "nondim_year"
        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        p = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        a = rootgrp.createVariable('a', 'f4', ('a',),zlib=True)
        x.units = "nondim_m"
        a.units = "nondim_strength"
        p.units = "nondim_mmtoyear"
        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        a[:] = a_range
        p[:] = p_range
        time[:] = t_range
        print "Setting up 1D variables"
        rootgrp.createVariable('b1', 'f4', ('p','a','time', 'x',),zlib=True)
        rootgrp.createVariable('b2', 'f4', ('p','a','time', 'x',),zlib=True)
        rootgrp.createVariable('w' , 'f4', ('p','a','time', 'x',),zlib=True)
    print "Output: netCDF file was created: ", fname+".nc"

def save_p_a_snapshot(fname,p_i,a_i,p,a,t,sol,Es):
    """ Save snapshot of the four fields b1,b2,w, together with the time
    """
    assert (len(Es['n']) == 1), "Able to save p to a scan only for 1D simulations"
    b1,b2,w = np.split(sol,3,axis=1)
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['b1'][p_i,a_i,:,:] = b1
        rootgrp['b2'][p_i,a_i,:,:] = b2
        rootgrp['w'][p_i,a_i,:,:]  = w

def load_p_a_info(fname):
    with netCDF4.Dataset("%s.nc"%fname, 'r') as rootgrp:
        p = rootgrp['p'][:]
        a = rootgrp['a'][:]
        time = rootgrp['time'][:]
        x  = rootgrp['x'][:]
    return p,a,time,x

def load_p_a_grid(fname,pstep,astep):
    """ Load results from p to a grid scan
    """
    with netCDF4.Dataset("%s.nc"%fname, 'r') as rootgrp:
        b1  = rootgrp['b1'][pstep,astep,:,:]
        b2  = rootgrp['b2'][pstep,astep,:,:]
        w   = rootgrp['w'][pstep,astep,:,:]
        a   = rootgrp['a'][astep]
        p   = rootgrp['p'][pstep]
    return b1,b2,w,p,a

def plot_p_a_grid(fname,pstep,astep,t=None):
    import matplotlib.pylab as plt
    b1,b2,w,p,a=load_p_a_grid(fname,pstep,astep)
    fig,ax = plt.subplots(1,2)
    if t is None:
        t = [0,1]
    ax[0].imshow(b1,extent=[0,1,t[-1],t[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
    ax[1].imshow(b2,extent=[0,1,t[-1],t[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
    plt.show()

def make_p_a_movie(fname,a,fps=None):
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    ax1.set_title(r'$b_1$')
    ax2.set_title(r'$b_2$')
    ax1.set_ylabel(r'$t$')
    ax1.set_xlabel(r'$x$')
    ax2.set_xlabel(r'$x$')
    # row and column sharing
    ims=[]
    p_arr,a_arr,time,x=load_p_a_info(fname)
    if fps is None:
        fps = int(len(p_arr)/3)
    a_i = (np.abs(a_arr - a)).argmin()
    for p_i,p in enumerate(p_arr):
        b1,b2,w,p,a=load_p_a_grid(fname,p_i,a_i)
        im1=ax1.imshow(b1,extent=[0,1,time[-1],time[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
        im2=ax2.imshow(b2,extent=[0,1,time[-1],time[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
        ims.append([im1,im2])
    ani = animation.ArtistAnimation(fig, ims, interval=50,
                                    blit=False, repeat_delay=1000)
    ani.save('%s.mp4'%fname, fps=fps)
    print "Movie made! File name is %s.mp4"%fname
    
def save_sim_snapshot(fname,step,time,Vs,Es):
    """ Save snapshot of the four fields b1,b2,b3,s1,s2, together with the time
    """
    b1,b2,b3,s1,s2 = Vs[0],Vs[1],Vs[2],Vs[3],Vs[4]
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['time'][step] = time
        rootgrp['b1_coverage'][step] = float(len(np.argwhere(b1>(np.amax(b1)/2.0))))/float(b1.size)
        rootgrp['b2_coverage'][step] = float(len(np.argwhere(b2>(np.amax(b2)/2.0))))/float(b2.size)
        rootgrp['b3_coverage'][step] = float(len(np.argwhere(b3>(np.amax(b3)/2.0))))/float(b3.size)
        if len(Es['n']) == 1:
            rootgrp['b1'][step,:] = b1
            rootgrp['b2'][step,:] = b2
            rootgrp['b3'][step,:] = b3
            rootgrp['s1'][step,:] = s1
            rootgrp['s2'][step,:] = s2
        elif len(Es['n']) == 2:
            rootgrp['b1'][step,:,:] = b1
            rootgrp['b2'][step,:,:] = b2
            rootgrp['b3'][step,:,:] = b3
            rootgrp['s1'][step,:,:] = s1
            rootgrp['s2'][step,:,:] = s2

def movie_maker(fname):
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
            lx = float(getattr(rootgrp['setup'],'lx'))
            fig, ax = plt.subplots(1, 1,figsize=(3*1.618,3))
    #        axes = plt.gca()
            ax.set_ylim([-0.1,1.0])
            ax.set_xlim([0,lx])
            ax.set_title(r'$b_1$ in green and $b_2$ in green')
            x  = rootgrp['x'][:]
            b1 = rootgrp['b1'][:,:]
            b2 = rootgrp['b2'][:,:]
    #        w =  rootgrp['w'][:,:]
            for i in xrange(len(t)):
                line1, = ax.plot(x,b1[i],'g-',label=r'$b_1$')
                line2, = ax.plot(x,b2[i],'r-',label=r'$b_2$')
                ims.append([line1,line2])
        elif nd == 2:
            fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
            fig.subplots_adjust(right=0.8)
            ax1.set_aspect('equal', 'datalim')
            ax1.set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax1.set_title(r'$b_1$', fontsize=14)
            ax2.set_aspect('equal', 'datalim')
            ax2.set_adjustable('box-forced')
    #        ax2.autoscale(False)
            ax2.set_title(r'$b_2$', fontsize=14)
            b1 = rootgrp['b1'][:,:,:]
            b2 = rootgrp['b2'][:,:,:]
    #        w  = rootgrp['w'][:,:,:]
            for i in xrange(len(t)):
                im1 = ax1.imshow(b1[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                im2 = ax2.imshow(b2[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                ims.append([im1,im2])
                cbar_ax2 = fig.add_axes([0.85, 0.54, 0.03, 0.35])
                fig.colorbar(im2, cax=cbar_ax2)
        ani = animation.ArtistAnimation(fig, ims, interval=50, blit=False, repeat_delay=1000)
        ani.save('%s.mp4'%fname)
        print "Movie made! File name is %s.mp4"%fname

def create_animation(fname,showtime=False):
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
            ax[0].set_ylabel(r'$b_1$', fontsize=25)
            ax[1].set_ylabel(r'$b_2$', fontsize=25)
            ax[2].set_ylabel(r'$b_3$', fontsize=25)
            ax[2].set_xlabel(r'$x$', fontsize=25)
            x  = rootgrp['x'][:]
            ax[0].set_xlim([x[0],x[-1]])
            b1 = rootgrp['b1'][:,:]
            b2 = rootgrp['b2'][:,:]
            b3 = rootgrp['b3'][:,:]
            for i in xrange(len(t)):
                line1, = ax[0].plot(x,b1[i],'g-')
                line2, = ax[1].plot(x,b2[i],'g-')
                line3, = ax[2].plot(x,b3[i],'g-')
                if showtime:
                    ax.set_title(r'$b$ at $t={:4.3f}$'.format(t[i]), fontsize=25)
                ims.append([line1,line2,line3])
        elif nd == 2:
            fig, ax = plt.subplots(1,2,sharex=True,sharey=True)
            fig.subplots_adjust(right=0.8)
            ax[0].set_aspect('equal', 'datalim')
            ax[0].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[0].set_title(r'$b_1$', fontsize=25)
            b1 = rootgrp['b1'][:,:,:]
            ax[1].set_aspect('equal', 'datalim')
            ax[1].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[1].set_title(r'$b_1$', fontsize=25)
            b2 = rootgrp['b2'][:,:,:]
            ax[2].set_aspect('equal', 'datalim')
            ax[2].set_adjustable('box-forced')
    #        ax1.autoscale(False)
            ax[2].set_title(r'$b_1$', fontsize=25)
            b3 = rootgrp['b3'][:,:,:]
            for i in xrange(len(t)):
                im1, = ax.imshow(b1[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                im2, = ax.imshow(b2[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                im3, = ax.imshow(b3[i],cmap=plt.cm.YlGn, animated=True,vmin=0.0,vmax=1.0)
                if showtime:
                    ax.set_title(r'$b$ at $t={:4.3f}$'.format(t[i]), fontsize=25)
                ims.append([im1,im2,im3])
                cbar_ax2 = fig.add_axes([0.85, 0.35, 0.05, 0.55])
                fig.colorbar(im1, cax=cbar_ax2)
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
    ani.save('%s.mp4'%fname)
    print "Movie made!"

def image_maker(fname,parameters,file_format="png",only_b=True,title=None):
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plt
    length=parameters['l']
    if title is None:
        title = r"$b2s2$ Model for $p=$"+str(parameters["p"])
    # row and column sharing
    rootgrp = netCDF4.Dataset("%s.nc"%fname, 'r', format='NETCDF4')
#    t = rootgrp['time'][:]
    if len(parameters['n']) == 1:
        x  = rootgrp['x'][:]
        b1 = rootgrp['b1'][:,:]
        b2 = rootgrp['b2'][:,:]
        s1 = rootgrp['s1'][:,:]
        s2 = rootgrp['s2'][:,:]
        rootgrp.close()
        if only_b:
            fig, (ax1,ax2) = plt.subplots(2,1, sharex='col', sharey='row')
            ax1.set_ylim([0.0,1.0])
            ax2.set_ylim([0.0,1.0])
            ax1.set_xlim([0.0,length[0]])
            ax2.set_xlim([0.0,length[0]])
            ax1.plot(x,b1[-1],'g-', fontsize=14)
            ax1.set_ylabel(r"$b_1$")
            ax2.plot(x,b2[-1],'g-')
            ax2.set_ylabel(r"$b_2$", fontsize=14)
        else:
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 2, sharex='col', sharey='row')
    #        axes = plt.gca()
            ax1[0].set_ylim([0.0,1.0])
            ax2[0].set_ylim([0.0,1.0])
            ax3[0].set_ylim([0.0,parameters['s_fc']])
            ax4[0].set_ylim([0.0,parameters['s_fc']])
            ax1[1].set_ylim([0.0,1.0])
            ax2[1].set_ylim([0.0,1.0])
            ax3[1].set_ylim([0.0,parameters['s_fc']])
            ax4[1].set_ylim([0.0,parameters['s_fc']])
            ax1[0].plot(x,b1[0],'g-')
            ax1[1].plot(x,b1[-1],'g-')
            ax2[0].plot(x,b2[0],'g-')
            ax2[1].plot(x,b2[-1],'g-')
            ax3[0].plot(x,s1[0],'b-')
            ax3[1].plot(x,s1[-1],'b-')
            ax4[0].plot(x,s2[0],'b-')
            ax4[1].plot(x,s2[-1],'b-')
    if len(parameters['n']) == 2:
        b1 = rootgrp['b1'][:,:,:]
        b2 = rootgrp['b2'][:,:,:]
        s1 = rootgrp['s1'][:,:,:]
        s2 = rootgrp['s2'][:,:,:]
        rootgrp.close()
        if only_b:
            fig, (ax1,ax2) = plt.subplots(1,2, sharex='col', sharey='row')
            fig.subplots_adjust(right=0.8)
            ax1.set_aspect('equal', 'datalim')
            ax1.set_adjustable('box-forced')
            ax1.set_title(r"$b_1$", fontsize=14)            
            ax2.set_aspect('equal', 'datalim')
            ax2.set_adjustable('box-forced')
            ax2.set_title(r"$b_2$", fontsize=14)            
            ax1.imshow(b1[-1],cmap=plt.cm.YlGn,vmin=0.0,vmax=1.0)
            im2=ax2.imshow(b2[-1],cmap=plt.cm.YlGn,vmin=0.0,vmax=1.0)
            cbar_ax2 = fig.add_axes([0.85, 0.3, 0.03, 0.35])
            fig.colorbar(im2, cax=cbar_ax2)
        else:
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 2, sharex='col', sharey='row')
            fig.subplots_adjust(right=0.8)
            ax1[0].set_aspect('equal', 'datalim')
            ax1[0].set_adjustable('box-forced')
            ax2[0].set_aspect('equal', 'datalim')
            ax2[0].set_adjustable('box-forced')
            ax3[0].set_aspect('equal', 'datalim')
            ax3[0].set_adjustable('box-forced')
            ax4[0].set_aspect('equal', 'datalim')
            ax4[0].set_adjustable('box-forced')
            ax1[1].set_aspect('equal', 'datalim')
            ax1[1].set_adjustable('box-forced')
            ax2[1].set_aspect('equal', 'datalim')
            ax2[1].set_adjustable('box-forced')
            ax3[1].set_aspect('equal', 'datalim')
            ax3[1].set_adjustable('box-forced')
            ax4[1].set_aspect('equal', 'datalim')
            ax4[1].set_adjustable('box-forced')
            im1a = ax1[0].imshow(b1[0],cmap=plt.cm.YlGn, vmin=0.0,vmax=1.0)
            im1b = ax1[1].imshow(b1[-1],cmap=plt.cm.YlGn,vmin=0.0,vmax=1.0)
            im2a = ax2[0].imshow(b2[0],cmap=plt.cm.YlGn, vmin=0.0,vmax=1.0)
            im2b = ax2[1].imshow(b2[-1],cmap=plt.cm.YlGn,vmin=0.0,vmax=1.0)
            im3a = ax3[0].imshow(s1[0],cmap=plt.cm.Blues, vmin=0.0,vmax=parameters['s_fc'])
            im3b = ax3[1].imshow(s1[-1],cmap=plt.cm.Blues,vmin=0.0,vmax=parameters['s_fc'])
            im4a = ax4[0].imshow(s2[0],cmap=plt.cm.Blues, vmin=0.0,vmax=parameters['s_fc'])
            im4b = ax4[1].imshow(s2[-1],cmap=plt.cm.Blues,vmin=0.0,vmax=parameters['s_fc'])
            cbar_ax2 = fig.add_axes([0.85, 0.54, 0.03, 0.35])
            fig.colorbar(im2a, cax=cbar_ax2)
            cbar_ax4 = fig.add_axes([0.85, 0.1, 0.03, 0.35])
            fig.colorbar(im3a, cax=cbar_ax4)
#    plt.title(title)
    plt.savefig('%s.png'%fname,format=file_format)
    print "Plots produced!"

def last_frame_arrays(fname,parameters,file_format="npy"):
    rootgrp = netCDF4.Dataset("%s.nc"%fname, 'r', format='NETCDF4')
    if len(parameters['n']) == 1:
        b1 = rootgrp['b1'][:,:]
        b2 = rootgrp['b2'][:,:]
        s1 = rootgrp['s1'][:,:]
        s2 = rootgrp['s2'][:,:]
    if len(parameters['n']) == 2:
        b1 = rootgrp['b1'][:,:,:]
        b2 = rootgrp['b2'][:,:,:]
        s1 = rootgrp['s1'][:,:,:]
        s2 = rootgrp['s2'][:,:,:]
    data = np.array([b1[-1],b2[-1],s1[-1],s2[-1]])
    if file_format == "npy":
        np.save(fname+".npy", data.T)
        print "Last frame saved in file ", fname+".npy"
    elif file_format == "txt":
        np.savetxt(fname+".dat", data.T, fmt='%10.8f', delimiter=' ', header='B1 B2 S1 S2')
        print "Last frame saved in file ", fname+".dat"
    elif file_format == "auto":
        x = np.linspace(0,1,len(b1[-1]))
        b1x = np.gradient(b1[-1])
        b2x = np.gradient(b2[-1])
        s1x = np.gradient(s1[-1])
        s2x = np.gradient(s2[-1])
        data=np.vstack((x,data,b1x,b2x,s1x,s2x))
        np.savetxt(fname+".dat", data.T, fmt='%10.8f', delimiter=' ')
        print "Last frame saved to AUTO format file ", fname+".dat"
    else:
        print "Uknown format, frame not saved"
    rootgrp.close()

def close_file(rootgrp):
    rootgrp.close()
    print "Output: netCDF4 file closed."
    
def setup_continuation(Ps,Es, fname,npoints):
    """
    Opening an netCDF4 file with 3 groups:
    1 parameters - contains model's parameters set
    2 snapshots - will contain subgroups named after the continuation parameter
    value
    3 continuation - will contain datasets of array of stats values along the
    """
    print "Configuring netCDF4 file."
    rootgrp = netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4')
    parms = rootgrp.createGroup('parameters')
    if len(Es['n'])==1:
        setattr(parms,'nd',1)
        nx=Es['n'][0]
        lx=Es['l'][0]
        setattr(parms,'nx',int(nx))
        setattr(parms,'lx',lx)
    elif len(Es['n'])==2:
        setattr(parms,'nd',2)
        nx,ny=Es['n']
        lx,ly=Es['l']
        setattr(parms,'nx',int(nx))
        setattr(parms,'ny',int(ny))
        setattr(parms,'lx',float(lx))
        setattr(parms,'ly',float(ly))
    for k,v in Es.items():
        if v==None or k=='l' or k=='n':
            pass
        elif type(v)!=str:
            setattr(parms, k, float(v))
        elif type(v)==str:
            setattr(parms, k, str(v))
    setattr(parms,'conv_P',Ps['conv_P'])
    setattr(parms,'K1',Ps['K1'])
    setattr(parms,'K2',Ps['K2'])
    rootgrp.description = "Continuation dataset for b2s1 model."
    rootgrp.history = "Created " + t.ctime(t.time())
    xdim = rootgrp.createDimension("x", Es['n'][0])
    prec_dim = rootgrp.createDimension('prec', npoints)
    prec = rootgrp.createVariable('prec', 'f8', ('prec',),zlib=True)
    prec.units = "mm/year"
    x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
    x.units = "m"
    x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
    if len(Es['n']) == 1:
        print "Setting up 1D variables"
        rootgrp.createVariable('b1', 'f8', ('prec', 'x',),zlib=True)
        rootgrp.createVariable('b2', 'f8', ('prec', 'x',),zlib=True)
        rootgrp.createVariable('s1', 'f8', ('prec', 'x', ),zlib=True)
        rootgrp.createVariable('s2', 'f8', ('prec', 'x', ),zlib=True)
    elif len(Es['n']) == 2:
        print "Setting up 2D variables"
        ydim = rootgrp.createDimension("y", Es['n'][1])
        y = rootgrp.createVariable('y', 'f4', ('y',),zlib=True)
        y.units = "m"
        y[:] = np.linspace(0,Es['l'][1], Es['n'][1])
        rootgrp.createVariable('b1', 'f8', ('prec', 'x', 'y',),zlib=True)
        rootgrp.createVariable('b2', 'f8', ('prec', 'x', 'y',),zlib=True)
        rootgrp.createVariable('s1', 'f8',  ('prec', 'x', 'y',),zlib=True)
        rootgrp.createVariable('s2', 'f8',  ('prec', 'x', 'y',),zlib=True)
#    fields = {'b1':b1,'b2':b2,'s':s,'h':h}
    print "Output: netCDF file was created: ", fname+".nc"
    return rootgrp

def save_continuation_snapshot(rootgrp, state, prec, noutput, parameters):
    """ Save snapshot of the four fields b1,b2,s,h, together with the time in
    years
    """
    b1,b2,s1,s2 = state[0],state[1],state[2],state[3]
    rootgrp['prec'][noutput] = prec
    if len(parameters['n']) == 1:
        rootgrp['b1'][noutput,:] = b1
        rootgrp['b2'][noutput,:] = b2
        rootgrp['s1'][noutput,:]  = s1
        rootgrp['s2'][noutput,:]  = s2
    elif len(parameters['n']) == 2:
        rootgrp['b1'][noutput,:,:] = b1
        rootgrp['b2'][noutput,:,:] = b2
        rootgrp['s1'][noutput,:,:]  = s1
        rootgrp['s2'][noutput,:,:]  = s2

def plot_sol_p(fname,Ps,Es,prec,file_format="png"):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # row and column sharing
    rootgrp = netCDF4.Dataset("%s.nc"%fname, 'r', format='NETCDF4')
    prec_array = rootgrp['prec'][:]
    indx = find_nearest(prec_array,prec)
    if len(Es['n']) == 1:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex='col', sharey='row')
#        axes = plt.gca()
        ax1[0].set_ylim([0.0,1.0])
        ax2[0].set_ylim([0.0,1.0])
        ax3[0].set_ylim([0.0,Ps['s_fc']])
        ax4[0].set_ylim([0.0,Ps['s_fc']])
        ax1[1].set_ylim([0.0,1.0])
        ax2[1].set_ylim([0.0,1.0])
        ax3[1].set_ylim([0.0,Ps['s_fc']])
        ax4[1].set_ylim([0.0,Ps['s_fc']])
        x  = rootgrp['x'][:]
        b1 = rootgrp['b1'][:,:]
        b2 = rootgrp['b2'][:,:]
        s1 = rootgrp['s1'][:,:]
        s2 = rootgrp['s2'][:,:]
        ax1.plot(x,b1[indx],'g-')
        ax2.plot(x,b2[indx],'g-')
        ax3.plot(x,s1[indx],'b-')
        ax4.plot(x,s2[indx],'b-')
    if len(Es['n']) == 2:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(2, 2, sharex='col', sharey='row')
        ax1.set_aspect('equal', 'datalim')
        ax1.set_adjustable('box-forced')
        ax2.set_aspect('equal', 'datalim')
        ax2.set_adjustable('box-forced')
        ax3.set_aspect('equal', 'datalim')
        ax3.set_adjustable('box-forced')
        ax4.set_aspect('equal', 'datalim')
        ax4.set_adjustable('box-forced')
        b1 = rootgrp['b1'][:,:,:]
        b2 = rootgrp['b2'][:,:,:]
        s1 = rootgrp['s1'][:,:,:]
        s2 = rootgrp['s2'][:,:,:]
        ax1.imshow(b1[indx],cmap=plt.cm.YlGn, vmin=0.0,vmax=1.0)
        ax2.imshow(b2[indx],cmap=plt.cm.YlGn, vmin=0.0,vmax=1.0)
        ax3.imshow(s1[indx],cmap=plt.cm.Blues, vmin=0.0,vmax=1.0)
        ax4.imshow(s2[indx],cmap=plt.cm.Blues, vmin=0.0,vmax=1.0)
#        plt.colorbar()
    plt.savefig('%s.png'%fname,format=file_format)
    print "Plots produced!"
    rootgrp.close()

def save_cont(fname):
    rootgrp = netCDF4.Dataset("%s.nc"%fname, 'r', format='NETCDF4')
    ndim=int(getattr(rootgrp['parameters'],'nd'))
    prec_array = rootgrp['prec'][:]
    if ndim == 1:
        b1sol = rootgrp['b1'][:,:]
        b2sol = rootgrp['b2'][:,:]
        s1sol = rootgrp['s1'][:,:]
        s2sol = rootgrp['s2'][:,:]
        b1=np.mean(b1sol,axis=-1)
        b2=np.mean(b2sol,axis=-1)
        s1=np.mean(s1sol,axis=-1)
        s2=np.mean(s2sol,axis=-1)        
    if ndim == 2:
        b1sol = rootgrp['b1'][:,:,:]
        b2sol = rootgrp['b2'][:,:,:]
        s1sol = rootgrp['s1'][:,:,:]
        s2sol = rootgrp['s2'][:,:,:]
        b1=np.mean(b1sol,axis=-1)
        b2=np.mean(b2sol,axis=-1)
        s1=np.mean(s1sol,axis=-1)
        s2=np.mean(s2sol,axis=-1)
        b1=np.mean(b1,axis=-1)
        b2=np.mean(b2,axis=-1)
        s1=np.mean(s1,axis=-1)
        s2=np.mean(s2,axis=-1)
    lx = getattr(rootgrp['parameters'],'lx')
    nx = getattr(rootgrp['parameters'],'nx')
    K1 = getattr(rootgrp['parameters'],'K1')
    K2 = getattr(rootgrp['parameters'],'K2')
    conv_P = getattr(rootgrp['parameters'],'conv_P')
    rootgrp.close()
    branch= {'p':prec_array,'b1sol':b1sol,'b2sol':b2sol,'s1sol':s1sol,'s2sol':s2sol,
            'b1':b1,'b2':b2,'s1':s1,'s2':s2,'lx':lx,
            'nx':nx,'nd':ndim,'K1':K1,'K2':K2,'conv_P':conv_P,'sol':True}
    import deepdish as dd
    dd.io.save(fname+'.hdf5', [branch])