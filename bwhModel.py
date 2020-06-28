# -*- coding: utf-8 -*-
"""
#  bwhModel3.py
#
#  Two layers model
#
#  Copyright 2016 Omer Tzuk <omertz@post.bgu.ac.il>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
"""

__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import time
from sympy import symbols, Matrix,lambdify
from sympy.utilities.autowrap import ufuncify
import numpy as np
import scipy.linalg as linalg
from scipy.integrate import odeint,solve_ivp
from scipy.fftpack import fftn, ifftn
import scipy.sparse as sparse
from utilities import handle_netcdf as hn
import deepdish.io as dd


parser = 'data/Parameter_sets/BWH_Ps_adaptive_dynamics_set13.csv'

Es_normal={'rhs':"BWH_adaptive_dynamics_v1",
        'n':(128,128),
        'l':(256.0,256.0),
        'bc':"neumann",
        'it':"pseudo_spectral",
        'dt':0.1,
        'analyze':False,
        'verbose':True,
        'setPDE':True}

def main():
    global m,p
    m = bwhModel(Es=Es_normal,Ps=parser,Vs="random")
    return 0

class bwhModel(object):
    def __init__(self,Ps,Es,Vs=None):
        if type(Ps)==str:
            self.Psfname=Ps
            if Ps.endswith('csv'):
                import pandas as pd
                self.p = pd.read_csv(Ps, index_col=0, squeeze=True).to_dict()
            elif Ps.endwith('hdf'):
                self.p=dd.load(Ps)
        else:
            self.Psfname=None
            self.p = Ps
        self.setup=Es
#        self.Vs=Vs
        self.verbose=Es['verbose']
        if self.verbose:
            start=time.time()
        self.set_equations()
        self.dt = 0.1
        self.time_elapsed = 0
        if self.setup['setPDE']:
            self.p['nd']=len(Es['n'])
            if self.p['nd']==2:
                self.p['nx'],self.p['ny']=Es['n']
                self.p['lx'],self.p['ly']=Es['l']
                self.l=[self.p['lx'],self.p['ly']]
                self.n=[self.p['nx'],self.p['ny']]
                self.dg  = tuple([l/float(n) for l,n in zip(self.l,self.n)])
                self.dx  = self.dg[0]
            elif self.p['nd']==1:
                self.dg=[Es['l'][0]/float(Es['n'][0])]
                self.dx=self.dg[0]
            self.dx2 = self.dx**2
            self.dt=Es['dt']*self.dx2 / np.amax(self.diffusion_coeffs)
            self.X = np.linspace(0,Es['l'][0],Es['n'][0])
            from utilities.laplacian_sparse import create_laplacian #,create_gradient
            self.lapmat=create_laplacian(self.setup['n'],self.setup['l'], self.setup['bc'] , self.diffusion_coeffs,verbose=self.verbose)
#            self.gradmat=create_gradient(self.setup['n'],self.setup['l'], self.setup['bc'] , [1.0,self.p['Dw'],self.p['Dh']])
            self.set_integrator()
            if self.verbose:
                print("Laplacian created")
        if Vs is not None:
            self.setup_initial_condition(Vs)
        if self.verbose:
            print("Time to setup: ",time.time()-start)
    """ Setting up model equations """
    def set_equations(self):
        t = symbols('t')
        self.Ps_symbols={}
        for key in list(self.p.keys()):
            self.Ps_symbols[key] = symbols(key)
        p=self.Ps_symbols
        if self.setup['rhs']=="BWH_adaptive_dynamics_v1":
            """ oz with tradeoff between E and K, with Q proportional to K
            Meaning that the higher K is, the Q is higher, and the Infitration
            is lower.
            """
            B,W,H = symbols('B W H')
            self.Vs_symbols = [B,W,H]
            self.setup['Vs_symbols'] = B,W,H
            self.setup['nvar']=len(self.setup['Vs_symbols'])
            self.diffusion_coeffs = [self.p['DB'],self.p['DW'],self.p['DH']]
            K=p['Kmax']+p['chi']*(p['Kmin']-p['Kmax'])
            M=p['Mmax']+p['chi']*(p['Mmin']-p['Mmax'])
            Y=p['Ymax']+p['chi']*(p['Ymin']-p['Ymax'])
            I = p['A']*((Y*B+p['f']*p['Q'])/(Y*B+p['Q']))
            self.dBdt_eq=p['Lambda']*(1.0-(B/(B+K)))*W*B-M*B
            self.dWdt_eq=I*H-((p['N']*W)/(1.0+p['R']*B))-p['Gamma']*W*B
            self.dHdt_eq=p['P']-I*H
            self.Sm_eq = p['Lambda']*(1.0-(B/(B+K)))*W-M
            self.Sm = ufuncify([B,W,H,t,p['P'],p['chi']],[self.sub_parms(self.Sm_eq)])
            self.dSmdchi_eq = self.Sm_eq.diff(p['chi'])
            self.dSmdchi = ufuncify([B,W,H,t,p['P'],p['chi']],[self.sub_parms(self.dSmdchi_eq)])
        """ Creating numpy functions """
        symeqs = Matrix([self.dBdt_eq,self.dWdt_eq,self.dHdt_eq])
        # For ODE computations
        self.ode  = lambdify((B,W,H,t,p['P'],p['chi']),self.sub_parms(symeqs),"numpy",dummify=False)
        # For PDE computations
        self.dBdt = ufuncify([B,W,H,t,p['P'],p['chi']],[self.sub_parms(self.dBdt_eq)])
        self.dWdt = ufuncify([B,W,H,t,p['P'],p['chi']],[self.sub_parms(self.dWdt_eq)])
        self.dHdt = ufuncify([B,W,H,t,p['P'],p['chi']],[self.sub_parms(self.dHdt_eq)])
        localJac   = symeqs.jacobian(Matrix([B,W,H]))
        self.sym_localJac = localJac
        self.localJac = lambdify((B,W,H,t,p['P'],p['chi']),self.sub_parms(localJac),"numpy",dummify=False)
        if self.setup['setPDE'] and self.setup['analyze']:
            self.dbdb = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[0,0])])
            self.dbdw = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[0,1])])
            self.dbdh = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[0,2])])
            self.dwdb = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[1,0])])
            self.dwdw = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[1,1])])
            self.dwdh = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[1,2])])
            self.dhdb = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[2,0])])
            self.dhdw = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[2,1])])
            self.dhdh = ufuncify([B,W,H,p['chi']],[self.sub_parms(localJac[2,2])])
        if self.verbose:
            self.print_equations()
            print ("Local Jacobian:" ,localJac)
            if self.setup['setPDE'] and self.setup['analyze']:
                print ("Linear analysis Jacobian: ", jaclinanalysis)

    """ Printing and parameters related functions """
    def print_parameters(self):
        print (self.p)
    def print_equations(self,numeric=False):
        if numeric:
            print ("dBdt = ", self.sub_parms(self.dBdt_eq))
            print ("dWdt = ", self.sub_parms(self.dWdt_eq))
            print ("dHdt = ", self.sub_parms(self.dHdt_eq))
        else:
            print ("dBdt = ", self.dBdt_eq)
            print ("dWdt = ", self.dWdt_eq)
            print ("dHdt = ", self.dHdt_eq)
    def print_latex_equations(self):
        from sympy import latex
        print ("\partial_t B = ",latex(self.dBdt_eq))
        print ("\partial_t W = ",latex(self.dWdt_eq))
        print ("\partial_t H = ",latex(self.dHdt_eq))
    """ Functions for use with scipy methods """
    def local(self,u ,t=0):
        B,W,H = u[0],u[1],u[2]
        return self.dBdt(B,W,H),self.dWdt(B,W,H),self.dHdt(B,W,H)
    def calc_ode_eigs(self,B,W,H,t=0,p=None,chi=None,beta=None,a=None):
        if p is not None:
            self.p['P']=p
        if chi is not None:
            self.p['chi']=chi
        if beta is not None:
            self.p['beta']=beta
        if a is not None:
            self.p['a']=a
        return linalg.eigvals(self.localJac(B,W,H,t,self.p['P'],self.p['chi']))
    def calc_SpatODE_eigs(self,B,W,H):
        return linalg.eigvals(self.SpatODEjac(B,W,H))

    def sigma_k_scan(self,B,W,H,k_range=[0,1.0],n=1000):
        k_range = np.linspace(k_range[0],k_range[1],n)
        MaxReSigma = np.zeros(n)
        MaxImSigma = np.zeros(n)
        for i,k in enumerate(k_range):
            eigsvalues=linalg.eigvals(self.jaclinanalysis(B,W,H,k))
            MaxReSigma[i]=np.amax(np.real(eigsvalues))
            MaxImSigma[i]=np.imag(eigsvalues[np.argmax(np.real(eigsvalues))])
        return np.array([k_range,MaxReSigma,MaxImSigma])

    def calc_linear_stability_analysis(self,B,W,H,k_range=[0,1.0],n=1000):
        k_scan = self.sigma_k_scan(B,W,H,k_range=[0,0.1],n=1000)
        return k_scan[0][np.argmax(k_scan[1])],np.amax(k_scan[1])
    """ Utilities """
    def sub_parms(self,eqs):
        B,W,H,t = symbols('B W H t')
        for key in list(self.p.keys()):
#            print key
            if key!='P' and key!='chi':
                eqs=eqs.subs(self.Ps_symbols[key],self.p[key])
        return eqs
    def calc_P_PET(self,p,w):
        P = self.p['dimpar']
        W = P['M0']/P['Lambda_max']*w
        PET = W*(P['NW']+P['Gamma']*(1+P['E']*P['K'])**2)
        PREC = p*self.p['conv_p_to_P']
        return PREC,PET


    """ Spatial functions """
    def set_integrator(self):
        integrator_type = {}
        integrator_type['scipy'] = self.scipy_integrate
        integrator_type['euler'] = self.euler_integrate
        integrator_type['rk4'] = self.rk4_integrate
        integrator_type['pseudo_spectral'] = self.pseudo_spectral_integrate
        try:
            self.integrator = integrator_type[self.setup['it']]
        except KeyError:
            raise  ValueError("No such integrator : %s."%self.setup['it'])
        if self.setup['it']=='pseudo_spectral':
            self.dt*=100.0

    def rhs_pde(self,state,t=0):
        B,W,H=np.split(state,self.setup['nvar'])
        return np.ravel((self.dBdt(B,W,H,t,self.p['P'],self.p['chi']),
                         self.dWdt(B,W,H,t,self.p['P'],self.p['chi']),
                         self.dHdt(B,W,H,t,self.p['P'],self.p['chi']))) + self.lapmat*state

    def rhs_pde_nonlindiff(self,state,t=0):
        B,W,H=np.split(state,self.setup['nvar'])
        hsq=h**2
        state_hsq = np.concatenate((B,W,Hsq))
        return np.ravel((self.dBdt(B,W,H,t,self.p['P'],self.p['chi']),
                         self.dWdt(B,W,H,t,self.p['P'],self.p['chi']),
                         self.dHdt(B,W,H,t,self.p['P'],self.p['chi']))) + self.lapmat*state_hsq

    def rhs_ode(self,state,t=0):
        B,W,H=state
        return self.ode(B,W,H,t,self.p['P'],self.p['chi']).T[0]
    def scipy_ode_rhs(self,t,state):
        B,W,H=state
        return np.squeeze(self.ode(B,W,H,t,self.p['P'],self.p['chi']))
    def scipy_ode_jac(self,t,state):
        B,W,H=state
        return self.localJac(B,W,H,t,self.p['P'],self.p['chi'])
    def calc_pde_analytic_jacobian(self,state):
        B,W,H=np.split(state,self.setup['nvar'])
        dbdb= sparse.diags(self.dbdb(B,W,H,self.p['chi']))
        dbdw= sparse.diags(self.dbdw(B,W,H,self.p['chi']))
        dbdh= sparse.diags(self.dbdh(B,W,H,self.p['chi']))
        dwdb= sparse.diags(self.dwdb(B,W,H,self.p['chi']))
        dwdw= sparse.diags(self.dwdw(B,W,H,self.p['chi']))
        dwdh= sparse.diags(self.dwdh(B,W,H,self.p['chi']))
        dhdb= sparse.diags(self.dhdb(B,W,H,self.p['chi']))
        dhdw= sparse.diags(self.dhdw(B,W,H,self.p['chi']))
        dhdh= sparse.diags(self.dhdh(B,W,H,self.p['chi']))
        local  = sparse.bmat([[dbdb,dbdw,dbdh],
                              [dwdb,dwdw,dwdh],
                              [dhdb,dhdw,dhdh]])
        return sparse.csc_matrix(local)+sparse.csc_matrix(self.lapmat)

    def calc_ode_numerical_jacobian(self,B,W,H,delta=0.00000001):
        state = np.array([B,W,H])
        jacobian = []
        for j in range(len(state)):
            state_plus = np.copy(state)
            state_minus = np.copy(state)
            state_plus[j] = state_plus[j]+delta
            state_minus[j] = state_minus[j]-delta
            jacobian.append((np.array(self.dudt(state_plus))-np.array(self.dudt(state_minus)))/(2.0*delta))
        return np.array(jacobian).T
    def check_pde_jacobians(self,n=100):
        import time
        timer_analytic=0
        timer_numeric=0
        error = np.zeros(n)
        for i in range(n):
            print (i)
            b=np.random.random(self.setup['n'])
            w=np.random.random(self.setup['n'])
            h=np.random.random(self.setup['n'])
            state=np.ravel((B,W,H))
            start_time=time.time()
            numeric=self.calc_pde_numerical_jacobian(state)
            mid_time=time.time()
            analytic=self.calc_pde_analytic_jacobian(state)
            end_time=time.time()
            timer_numeric+=(mid_time-start_time)
            timer_analytic+=(end_time-mid_time)
            error[i]=np.max(np.abs(numeric-analytic))
        print ("Max difference is ",np.max(error), ", and mean difference is ",np.mean(error))
        print ("Average speed for numeric ", timer_numeric/float(n))
        print ("Average speed for analytic ", timer_analytic/float(n))
        print ("Analytic ", float(timer_numeric)/float(timer_analytic)," faster.")

    def calc_pde_numerical_jacobian(self,state,delta=0.00000001):
        n = len(state)
        jacobian = []
        for j in range(n):
            state_plus = np.copy(state)
            state_minus = np.copy(state)
            state_plus[j] = state_plus[j]+delta
            state_minus[j] = state_minus[j]-delta
            jacobian.append((self.rhs_pde(state_plus)-self.rhs_pde(state_minus))/(2.0*delta))
        return np.array(jacobian).T

    def calc_numeric_pde_eigs(self,state):
        return linalg.eigvals(self.calc_pde_numerical_jacobian(state))
    def calc_analytic_pde_eigs(self,state):
        return sparse.linalg.eigs(self.calc_pde_analytic_jacobian(state),k=3)[0]

    def check_convergence(self,state,previous_state,tolerance=1.0e-5):
        return np.max(np.abs(state-previous_state))<tolerance

    def update_parameters(self,parameters):
        intersection=[i for i in list(self.p.keys()) if i in parameters]
        if intersection:
            for key in intersection:
                if self.setup['verbose']:
                    print (str(key)+"="+str(parameters[key]))
                self.p[key]=parameters[key]

    def integrate(self,initial_state=None,step=10,
                  max_time = 1000,tol=1.0e-5,plot=False,savefile=None,
                  create_movie=False,check_convergence=True,
                  sim_step=None,**kwargs):
        if kwargs:
            self.update_parameters(kwargs)
        self.filename = savefile
        if initial_state is None:
            initial_state = self.initial_state
        self.time_elapsed=0
        if sim_step is None:
            self.sim_step=0
        else:
            self.sim_step=sim_step
        if savefile is not None:
            hn.setup_simulation(savefile,self.p,self.setup)
            hn.save_sim_snapshot(savefile,self.sim_step,self.time_elapsed,
                                 self.split_state(initial_state),self.setup)
#        old_result = initial_state
        converged=False
        result = []
        result.append(initial_state)
#        t,result = self.integrator(initial_state,p=p,step=10,finish=10,savefile=self.filename)
        if self.setup['verbose']:
            start=time.time()
            print ("Step {:4d}, Time = {:5.1f}".format(self.sim_step,self.time_elapsed))
        while not converged and self.time_elapsed<=max_time:
            old_result = result[-1]
            t,result = self.integrator(result[-1],step=step,finish=step)
#            self.time_elapsed+=t[-1]
            self.sim_step=self.sim_step+1
            if savefile is not None:
                hn.save_sim_snapshot(savefile,self.sim_step,self.time_elapsed,
                                     self.split_state(result[-1]),self.setup)
            if self.setup['verbose']:
                print ("Step {:4d}, Time = {:10.6f}, diff = {:7f}".format(self.sim_step,self.time_elapsed,np.max(np.abs(result[-1]-old_result))))
            if check_convergence:
                converged=self.check_convergence(result[-1],old_result,tol)
                if converged:
                    print ("Convergence detected")
        if self.setup['verbose']:
            print ("Integration time was {} s".format(time.time()-start))
        if savefile is not None and create_movie:
            print ("Creating movie...")
            if type(create_movie)==str:
                hn.create_animation(savefile,create_movie)
            else:
                hn.create_animation(savefile)
        return result[-1]

    """ Integrators step functions """
    def update_parameters(self,parameters):
        intersection=[i for i in self.p.keys() if i in parameters]
        if intersection:
            for key in intersection:
                if self.setup['verbose']:
                    print(str(key)+"="+str(parameters[key]))
                self.p[key]=parameters[key]

    def scipy_integrate(self,initial_state,step=0.1,finish=1000):
        """ """
#        print "Integration using scipy odeint"
        t = np.arange(0,finish+step, step)
        self.time_elapsed+=finish
        return t,odeint(self.rhs_pde, initial_state, t)
    def ode_integrate(self,initial_state,step=0.1,start=0,finish=1000,
                          method='BDF',**kwargs):
        """ Using the new scipy interface to BDF method for stiff equations
        with option to switch to other methods
        """
        if kwargs:
            self.update_parameters(kwargs)
        t = np.arange(start,finish+step, step)
        sol=solve_ivp(fun=self.scipy_ode_rhs,t_span=(start,finish),
                      y0=initial_state,method=method,
                      t_eval=t,jac=self.scipy_ode_jac)
        return sol.t,sol.y
    def ode_integrate_old(self,initial_state,step=0.1,start=0,finish=1000,**kwargs):
        """ """
        if kwargs:
            self.update_parameters(kwargs)
        t = np.arange(start,finish+step, step)
        return t,odeint(self.rhs_ode,initial_state, t).T

    def euler_integrate(self,initial_state,step=0.1,finish=1000):
        """ """
        time = np.arange(0,finish+step,step)
        result=[]
        t=0
        result.append(initial_state)
        for tout in time[1:]:
            old=result[-1]
            while t < tout:
#                print "t=",t
                new=old+self.dt*self.rhs_pde(old)
                old=new
                t+=self.dt
            result.append(old)
        return t,result
    def rk4_integrate(self,initial_state,step=0.1,finish=1000):
        """ """
#        print "Integration using rk4 step"
        time = np.arange(0,finish+step,step)
        result=[]
        t=0
        result.append(initial_state)
        for tout in time[1:]:
            old=result[-1]
            while t < tout:
                k1=self.dt*self.rhs_pde(old,self.time_elapsed)
                k2=self.dt*self.rhs_pde(old+0.5*k1,self.time_elapsed)
                k3=self.dt*self.rhs_pde(old+0.5*k2,self.time_elapsed)
                k4=self.dt*self.rhs_pde(old+k3,self.time_elapsed)
                new=old+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
                old=new
                t+=self.dt
                self.time_elapsed+=self.dt
            result.append(old)
        return time,result

    def pseudo_spectral_integrate(self,initial_state,step=0.1,finish=1000,**kwargs):
#        print "Integration using pseudo-spectral step"
        if kwargs:
            self.update_parameters(kwargs)
        time = np.arange(0,finish+step,step)
        result=np.zeros((len(time),len(initial_state)))
        t=0
        result[0]=initial_state
        for i,tout in enumerate(time[1:]):
            self.state=result[-1]
            B,W,H=self.state.reshape(self.setup['nvar'],*self.setup['n'])
            self.fftb=fftn(B)
            self.fftw=fftn(W)
            self.ffth=fftn(H)
            while t < tout:
                self.fftb = self.multb*(self.fftb + self.dt*fftn(self.dBdt(B,W,H,t,self.p['P'],self.p['chi'])))#.real
                self.fftw = self.multw*(self.fftw + self.dt*fftn(self.dWdt(B,W,H,t,self.p['P'],self.p['chi'])))#.real
                self.ffth = self.multh*(self.ffth + self.dt*fftn(self.dHdt(B,W,H,t,self.p['P'],self.p['chi'])))#.real
                B= ifftn(self.fftb).real
                W= ifftn(self.fftw).real
                H= ifftn(self.ffth).real
                t+=self.dt
                self.time_elapsed+=self.dt
            result[i+1]=np.ravel((B,W,H))
        self.state=result[-1]
        return time,result

    def spectral_multiplier(self,dt):
        n=self.setup['n']
        nx=n[0]
        dx=self.dx
        # wave numbers
        k=2.0*np.pi*np.fft.fftfreq(nx,dx)
        if len(n)==1:
            k2=k**2
        if len(n)==2:
            k2=np.outer(k,np.ones(nx))**2
            k2+=np.outer(np.ones(nx),k)**2
        # multiplier
        self.multb = np.exp(-dt*self.p['DB']*k2)
        self.multw = np.exp(-dt*self.p['DW']*k2)
        self.multh = np.exp(-dt*self.p['DH']*k2)

    def setup_initial_condition(self,Vs,**kwargs):
        n = self.setup['n']
        if type(Vs)==str:
            if Vs == "random":
                B = np.random.random(n)*0.5 + 0.1
                W= np.random.random(n)*(0.1) + 0.05
                H= np.random.random(n)*(0.1) + 0.05
            elif Vs == "bare":
                B = np.zeros(n)
                W= np.random.random(n)*(0.1) + 0.05
                H= np.random.random(n)*(0.1) + 0.05
            elif Vs == "uniform":
                p   = kwargs.get('p', None)
                chi = kwargs.get('chi', None)
                t,sol = self.ode_integrate([0.2,0.2,0.2],p=p,chi=chi)
                B= np.ones(n)*sol[0][-1]
                W= np.ones(n)*sol[1][-1]
                H= np.ones(n)*sol[2][-1]
            elif Vs == "tile":
                fields = kwargs.get('fields', None)
                B,W,H = np.split(fields,self.setup['nvar'])
                B  = np.tile(b,(self.setup['n'][0],1))
                W = np.tile(w,(self.setup['n'][0],1))
                H = np.tile(h,(self.setup['n'][0],1))
            self.initial_state = np.ravel((B,W,H))
        else:
            self.initial_state = Vs
        self.state = self.initial_state
        if self.setup['it'] == 'pseudo_spectral' and self.setup['setPDE']:
            self.spectral_multiplier(self.dt)
    """ Plot functions """

    def plotODEInt(self,p,a,initial_state=[0.8,0.2,0.2],
                   chi=[0,0.5,1.0],
                   step=0.1,start=0,finish=200):
        import matplotlib.pylab as plt
        plt.figure()
        cl = plt.cm.winter(np.linspace(0, 1, 3))
        chi = np.array(chi)
        for i,chi_v in enumerate(chi):
            t,sol = self.ode_integrate(initial_state,step,start,finish,p=p,a=a,chi=chi_v)
            plt.plot(t/self.p['conv_T_to_t'], sol[0],color=cl[i], lw=2,label=r'$\chi={:3.2f}$'.format(chi_v))
        plt.xlabel(r'$T \; [yr]$')
        plt.ylabel(r'$b$')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.show()


    def split_state(self,state):
        return state.reshape(self.setup['nvar'],*self.setup['n'])

    def plot(self,state=None,fontsize=12,update=False):
        if state is None:
            state=self.state
        import matplotlib.pylab as plt
        if update:
            plt.ion()
        B,W,H=state.reshape(self.setup['nvar'],*self.setup['n'])
        if len(self.setup['n'])==1:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    #        x = np.linspace(0,self.setup['l'][0],self.setup['n'][0])
            ax1.set_ylim([-0.1,1.0])
#            ax2.set_ylim([0.0,self.p['s_fc']])
#            ax3.set_ylim([0.0,self.p['s_fc']])
            ax1.plot(self.X,B)
            ax1.set_xlim([0,self.X[-1]])
            ax2.plot(self.X,W)
            ax3.plot(self.X,H)
            ax3.set_xlabel(r'$x$', fontsize=fontsize)
            ax1.set_ylabel(r'$B$', fontsize=fontsize)
            ax2.set_ylabel(r'$W$', fontsize=fontsize)
            ax3.set_ylabel(r'$H$', fontsize=fontsize)
        elif len(self.setup['n'])==2:
            fig, (ax1, ax2,ax3) = plt.subplots(1, 3, sharex=True, sharey=True)
            fig.subplots_adjust(right=0.8)
#            ax1.imshow(b,cmap=plt.cm.YlGn, vmin = bmin, vmax = bmax)
            ax1.imshow(B,cmap=plt.cm.YlGn)
            ax1.set_title(r'$B$', fontsize=fontsize)
            ax2.imshow(W,cmap=plt.cm.Blues)
            ax2.set_title(r'$W$', fontsize=fontsize)
#            ax3.imshow(w,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
            ax3.imshow(H,cmap=plt.cm.Blues)
            ax3.set_title(r'$H$', fontsize=fontsize)
#            im4=ax4.imshow(h,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
#            cbar_ax2 = fig.add_axes([0.85, 0.54, 0.03, 0.35])
#            fig.colorbar(im1, cax=cbar_ax2)
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()