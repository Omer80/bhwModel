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
#import netCDF4 as netcdf
from utilities import handle_netcdf3 as hn
#from scipy.io import netcdf
import deepdish.io as dd
#from tlm_parameters import par

Es_normal={'rhs':"oz_EQK",
        'n':(512,512),
        'l':(256.0,256.0),
        'bc':"neumann",
        'it':"pseudo_spectral",
        'dt':0.1,
        'analyze':False,
        'verbose':True,
        'setPDE':True}

def main():
    global m,p
    m = bwhModel(Es=Es_normal,Ps='auto/bwh_set4.hdf5',Vs="random")
    return 0

class bwhModel(object):
    def __init__(self,Ps,Es,Vs=None):
        if type(Ps)==str:
            self.Psfname=Ps
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
            self.dt=Es['dt']*self.dx2 / self.p['dh']
            self.X = np.linspace(0,Es['l'][0],Es['n'][0])
            from utilities.laplacian_sparse3 import create_laplacian #,create_gradient
            self.lapmat=create_laplacian(self.setup['n'],self.setup['l'], self.setup['bc'] , [1.0,self.p['dw'],self.p['dh']],verbose=self.verbose)
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
        if self.setup['rhs']=="oz":
            """ Ms version tradeoff of Lambda with Ms AND Gamma
            """
            b,w,h = symbols('b w h')
            self.setup['Vs_symbols'] = b,w,h
            self.setup['nvar']=len(self.setup['Vs_symbols'])
            from sympy.functions import cos as symcos
            g    = (1.0+p['eta']*b)*(1.0+p['eta']*b)
            evapw = (p['nuw']*w)/(1.0+p['rhow']*b)
            evaph = (p['nuh']*w)/(1.0+p['rhoh']*b)
            tras = p['gamma']*g*w*b
            i    = p['alpha']*((b+p['q']*p['f'])/(b+p['q']))
            self.p_t_eq = p['p']*(1.0+p['a']*symcos(2.0*np.pi*p['omegaf']*t/p['conv_T_to_t']))
            self.dbdt_eq = p['lamb_max']*g*w*b*(1.0-b)-b
            self.dwdt_eq = i*h-evapw-tras
            self.dhdt_eq = self.p_t_eq-evaph-i*h
        elif self.setup['rhs']=="tradeoff":
            """ Ms version tradeoff of Lambda with Ms AND Gamma 
            """
            b,w,h = symbols('b w h')
            self.Vs_symbols = b,w,h
            self.setup['nvar']=len(self.Vs_symbols)
            from sympy.functions import exp as sympexp
            from sympy.functions import cos as symcos
            sigma = 100
            g    = (1.0+p['eta']*b)*(1.0+p['eta']*b)
            evap = (p['nu']*w)/(1.0+p['rho']*b)
            midlamb = (p['lamb_max']+p['lamb_min'])/2.0
            self.lamb = p['lamb_max'] + (p['chi']**p['beta']) * (p['lamb_min'] - p['lamb_max'])
            self.gamma= p['gamma']*(self.lamb/midlamb)
            self.mu_s = p['mu_s_max'] + ((1.0-p['chi'])**p['beta']) * (0.0 - p['mu_s_max'])
            self.mu  = 1.0-self.mu_s*(1.0/(1.0 + sympexp(sigma*(h-(p['s_wp']+p['s_fos'])/2.0))))
#            ms   = 0.1+0.9*(1.0/(1.0 + sympexp(-sigma*(h-s_wp))))
            tras = self.gamma*g*w*b
            i    = p['alpha']*((b+p['q']*p['f'])/(b+p['q']))
            self.p_t_eq = p['p']*(1.0+p['a']*symcos(2.0*np.pi*p['omegaf']*t/p['conv_T_to_t']))
            self.dbdt_eq  = self.lamb*g*w*b*(1.0-b)-self.mu*b
            self.dwdt_eq = i*h-evap-tras
            self.dhdt_eq = self.p_t_eq-i*h
        elif self.setup['rhs']=="oz_EQK":
            """ oz with tradeoff between E and K, with Q proportional to K
            Meaning that the higher K is, the Q is higher, and the Infitration 
            is lower. 
            """
            b,w,h = symbols('b w h')
            self.Vs_symbols = [b,w,h]
            self.setup['Vs_symbols'] = b,w,h
            self.setup['nvar']=len(self.setup['Vs_symbols'])
            from sympy.functions import cos as symcos
#            from sympy.functions import exp as symexp
            q_min   = p['q']*(1.0-p['del_to'])
            q_max   = p['q']*(1.0+p['del_to'])
            eta_min = p['eta']*(1.0-p['del_to'])
            eta_max = p['eta']*(1.0+p['del_to'])
            K_min   = (1.0-p['del_to'])
            K_max   = (1.0+p['del_to'])
            K_to   = K_max    + (p['chi']**p['beta'])*(K_min-K_max)
            q_to   = (q_max   + (p['chi']**p['beta'])*(q_min-q_max))/K_to
            eta_to = (eta_max + ((1.0-p['chi'])**p['beta'])*(eta_min-eta_max))*K_to
            gam    = p['gamma']*K_to
            G = w*(1.0 + eta_to*b)*(1.0 + eta_to*b)
            I = (p['alpha']*((b + q_to*p['f'])/(b + q_to)))
            evapw = p['n']*((p['nuw'])/(1.0 + p['rhow']*b))*w
            evaph = p['n']*((p['nuh'])/(1.0 + p['rhoh']*b))*h
            tras  = gam*b*G
            self.p_t_eq = p['p']*(1.0+p['a']*symcos(2.0*np.pi*p['omegaf']*t/p['conv_T_to_t']))
            self.dbdt_eq=G*b*(1.0-b) - b
            self.dwdt_eq=I*h - evapw - tras
            self.dhdt_eq=self.p_t_eq - I*h - evaph
            self.conv_K = lambdify((b,p['chi']),self.sub_parms(b*K_to),"numpy",dummify=False)
        elif self.setup['rhs']=="oz_EQK_relax":
            """ Ms version tradeoff of Lambda with Ms AND Gamma
            """
            b,w,h = symbols('b w h')
            self.Vs_symbols = [b,w,h]
            self.setup['Vs_symbols'] = b,w,h
            self.setup['nvar']=len(self.setup['Vs_symbols'])
            from sympy.functions import cos as symcos
#            from sympy.functions import exp as symexp
            q_min   = p['q']*(1.0-p['del_to'])
            q_max   = p['q']*(1.0+p['del_to'])
            eta_min = p['eta']*(1.0-p['del_to'])
            eta_max = p['eta']*(1.0+p['del_to'])
            K_min   = (1.0-p['del_to'])
            K_max   = (1.0+p['del_to'])
            K_to    = K_max    + p['chi']*(K_min-K_max)
            q_to    = (q_max   + p['chi']*(q_min-q_max))/K_to
            eta_to  = (eta_max + (1.0-p['chi'])*(eta_min-eta_max))*K_to
            gam     = p['gamma']*K_to
            G = w*(1.0 + eta_to*b)*(1.0 + eta_to*b)
            I = (p['alpha']*((b + q_to*p['f'])/(b + q_to)))
            evapw = ((p['nuw'])/(1.0 + p['rhow']*b))*w
            evaph = ((p['nuh'])/(1.0 + p['rhoh']*b))*h
            tras  = gam*b*G
            self.dbdt_eq=G*b*(1.0-b) - b
            self.dwdt_eq=I*h - evapw - tras
            self.dhdt_eq= - I*h - evaph
            self.conv_K = lambdify((b,p['chi']),self.sub_parms(b*K_to),"numpy",dummify=False)
        """ Creating numpy functions """
        symeqs = Matrix([self.dbdt_eq,self.dwdt_eq,self.dhdt_eq])
        self.ode  = lambdify((b,w,h,t,p['p'],p['chi'],p['beta'],p['a'],p['omegaf']),self.sub_parms(symeqs),"numpy",dummify=False)
        self.dbdt = ufuncify([b,w,h,t,p['p'],p['chi'],p['beta'],p['a'],p['omegaf']],[self.sub_parms(self.dbdt_eq)])
        self.dwdt = ufuncify([b,w,h,t,p['p'],p['chi'],p['beta'],p['a'],p['omegaf']],[self.sub_parms(self.dwdt_eq)])
        self.dhdt = ufuncify([b,w,h,t,p['p'],p['chi'],p['beta'],p['a'],p['omegaf']],[self.sub_parms(self.dhdt_eq)])
        localJac   = symeqs.jacobian(Matrix([b,w,h]))
        self.sym_localJac = localJac
        self.localJac = lambdify((b,w,h,t,p['p'],p['chi'],p['beta'],p['a'],p['omegaf']),self.sub_parms(localJac),"numpy",dummify=False)
        if self.setup['setPDE'] and self.setup['analyze']:
            self.dbdb = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[0,0])])
            self.dbdw = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[0,1])])
            self.dbdh = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[0,2])])
            self.dwdb = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[1,0])])
            self.dwdw = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[1,1])])
            self.dwdh = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[1,2])])
            self.dhdb = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[2,0])])
            self.dhdw = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[2,1])])
            self.dhdh = ufuncify([b,w,h,p['chi']],[self.sub_parms(localJac[2,2])])
            k  = symbols('k')
            dw = symbols('dw')
            dh = symbols('dh')
            symeqs_lin_analysis = Matrix([self.dbdt_eq-b*k*k,self.dwdt_eq-w*dw*k*k,self.dhdt_eq-h*dh*k*k])
            jaclinanalysis = symeqs_lin_analysis.jacobian(Matrix([b,w,h]))
            self.symbolic_jaclinanalysis = jaclinanalysis
            self.jaclinanalysis = lambdify((b,w,h,k),self.sub_parms(jaclinanalysis),"numpy",dummify=False)
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
            print ("dbdt = ", self.sub_parms(self.dbdt_eq))
            print ("dwdt = ", self.sub_parms(self.dwdt_eq))
            print ("dhdt = ", self.sub_parms(self.dhdt_eq))
        else:
            print ("dbdt = ", self.dbdt_eq)
            print ("dwdt = ", self.dwdt_eq)
            print ("dhdt = ", self.dhdt_eq)
    def print_latex_equations(self):
        from sympy import latex
        print ("\partial_t b = ",latex(self.dbdt_eq))
        print ("\partial_t w = ",latex(self.dwdt_eq))
        print ("\partial_t h = ",latex(self.dhdt_eq))
    """ Functions for use with scipy methods """
    def local(self,u ,t=0):
        b,w,h = u[0],u[1],u[2]
        return self.dbdt(b,w,h),self.dwdt(b,w,h),self.dhdt(b,w,h)
    def calc_ode_eigs(self,b,w,h,t=0,p=None,chi=None,beta=None,a=None):
        if p is not None:
            self.p['p']=p
        if chi is not None:
            self.p['chi']=chi
        if beta is not None:
            self.p['beta']=beta
        if a is not None:
            self.p['a']=a
        return linalg.eigvals(self.localJac(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']))
    def calc_SpatODE_eigs(self,b,w,h):
        return linalg.eigvals(self.SpatODEjac(b,w,h))

    def sigma_k_scan(self,b,w,h,k_range=[0,1.0],n=1000):
        k_range = np.linspace(k_range[0],k_range[1],n)
        MaxReSigma = np.zeros(n)
        MaxImSigma = np.zeros(n)
        for i,k in enumerate(k_range):
            eigsvalues=linalg.eigvals(self.jaclinanalysis(b,w,h,k))
            MaxReSigma[i]=np.amax(np.real(eigsvalues))
            MaxImSigma[i]=np.imag(eigsvalues[np.argmax(np.real(eigsvalues))])
        return np.array([k_range,MaxReSigma,MaxImSigma])

    def calc_linear_stability_analysis(self,b,w,h,k_range=[0,1.0],n=1000):
        k_scan = self.sigma_k_scan(b,w,h,k_range=[0,0.1],n=1000)
        return k_scan[0][np.argmax(k_scan[1])],np.amax(k_scan[1])
    """ Utilities """
    def sub_parms(self,eqs):
        b,w,h,t = symbols('b w h t')
        for key in list(self.p.keys()):
#            print key
            if key!='p' and key!='chi' and key!='beta' and key!='a' and key!='omegaf':
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
        b,w,h=np.split(state,self.setup['nvar'])
        return np.ravel((self.dbdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']),
                         self.dwdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']),
                         self.dhdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']))) + self.lapmat*state
    
    def rhs_pde_nonlindiff(self,state,t=0):
        b,w,h=np.split(state,self.setup['nvar'])
        hsq=h**2
        state_hsq = np.concatenate((b,w,hsq))
        return np.ravel((self.dbdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']),
                         self.dwdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']),
                         self.dhdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']))) + self.lapmat*state_hsq

    def rhs_ode(self,state,t=0):
        b,w,h=state
        return self.ode(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']).T[0]
    def scipy_ode_rhs(self,t,state):
        b,w,h=state
        return np.squeeze(self.ode(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']))
    def scipy_ode_jac(self,t,state):
        b,w,h=state
        return self.localJac(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf'])
    def calc_pde_analytic_jacobian(self,state):
        b,w,h=np.split(state,self.setup['nvar'])
        dbdb= sparse.diags(self.dbdb(b,w,h,self.p['chi']))
        dbdw= sparse.diags(self.dbdw(b,w,h,self.p['chi']))
        dbdh= sparse.diags(self.dbdh(b,w,h,self.p['chi']))
        dwdb= sparse.diags(self.dwdb(b,w,h,self.p['chi']))
        dwdw= sparse.diags(self.dwdw(b,w,h,self.p['chi']))
        dwdh= sparse.diags(self.dwdh(b,w,h,self.p['chi']))
        dhdb= sparse.diags(self.dhdb(b,w,h,self.p['chi']))
        dhdw= sparse.diags(self.dhdw(b,w,h,self.p['chi']))
        dhdh= sparse.diags(self.dhdh(b,w,h,self.p['chi']))
        local  = sparse.bmat([[dbdb,dbdw,dbdh],
                              [dwdb,dwdw,dwdh],
                              [dhdb,dhdw,dhdh]])
        return sparse.csc_matrix(local)+sparse.csc_matrix(self.lapmat)

    def calc_ode_numerical_jacobian(self,b,w,h,delta=0.00000001):
        state = np.array([b,w,h])
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
            state=np.ravel((b,w,h))
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
    def pseudo_spectral_integrate(self,initial_state,step=0.1,finish=1000):
#        print "Integration using pseudo-spectral step"
        time = np.arange(0,finish+step,step)
        result=[]
        t=0
        result.append(initial_state)
        for tout in time[1:]:
            self.state=result[-1]
            b,w,h=self.state.reshape(self.setup['nvar'],*self.setup['n'])
            self.fftb=fftn(b)
            self.fftw=fftn(w)
            self.ffth=fftn(h)
            while t < tout:
                self.fftb = self.multb*(self.fftb + self.dt*fftn(self.dbdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf'])))#.real
                self.fftw = self.multw*(self.fftw + self.dt*fftn(self.dwdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf'])))#.real
                self.ffth = self.multh*(self.ffth + self.dt*fftn(self.dhdt(b,w,h,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf'])))#.real
                b= ifftn(self.fftb).real
                w= ifftn(self.fftw).real
                h= ifftn(self.ffth).real
                t+=self.dt
                self.time_elapsed+=self.dt
            self.state=np.ravel((b,w,h))
            self.sim_step+=1
            result.append(self.state)
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
        self.multb = np.exp(-dt*k2)
        self.multw = np.exp(-dt*self.p['dw']*k2)
        self.multh = np.exp(-dt*self.p['dh']*k2)

    def setup_initial_condition(self,Vs,**kwargs):
        n = self.setup['n']
        if type(Vs)==str:
            if Vs == "random":
                b = np.random.random(n)*0.5 + 0.1
                w= np.random.random(n)*(0.1) + 0.05
                h= np.random.random(n)*(0.1) + 0.05
            elif Vs == "bare":
                b = np.zeros(n)
                w= np.random.random(n)*(0.1) + 0.05
                h= np.random.random(n)*(0.1) + 0.05
            elif Vs == "uniform":
                p   = kwargs.get('p', None)
                chi = kwargs.get('chi', None)
                t,sol = self.ode_integrate([0.2,0.2,0.2],p=p,chi=chi)
                b= np.ones(n)*sol[0][-1]
                w= np.ones(n)*sol[1][-1]
                h= np.ones(n)*sol[2][-1]
            elif Vs == "tile":
                fields = kwargs.get('fields', None)
                b,w,h = np.split(fields,self.setup['nvar'])
                b  = np.tile(b,(self.setup['n'][0],1))
                w = np.tile(w,(self.setup['n'][0],1))
                h = np.tile(h,(self.setup['n'][0],1))                
            self.initial_state = np.ravel((b,w,h))
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
        b,w,h=state.reshape(self.setup['nvar'],*self.setup['n'])
        if len(self.setup['n'])==1:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    #        x = np.linspace(0,self.setup['l'][0],self.setup['n'][0])
            ax1.set_ylim([-0.1,1.0])
#            ax2.set_ylim([0.0,self.p['s_fc']])
#            ax3.set_ylim([0.0,self.p['s_fc']])
            ax1.plot(self.X,b)
            ax1.set_xlim([0,self.X[-1]])
            ax2.plot(self.X,w)
            ax3.plot(self.X,h)
            ax3.set_xlabel(r'$x$', fontsize=fontsize)
            ax1.set_ylabel(r'$b$', fontsize=fontsize)
            ax2.set_ylabel(r'$w$', fontsize=fontsize)
            ax3.set_ylabel(r'$h$', fontsize=fontsize)
        elif len(self.setup['n'])==2:
            fig, (ax1, ax2,ax3) = plt.subplots(1, 3, sharex=True, sharey=True)
            fig.subplots_adjust(right=0.8)
#            ax1.imshow(b,cmap=plt.cm.YlGn, vmin = bmin, vmax = bmax)
            ax1.imshow(b,cmap=plt.cm.YlGn)
            ax1.set_title(r'$b$', fontsize=fontsize)
            ax2.imshow(w,cmap=plt.cm.Blues)
            ax2.set_title(r'$h$', fontsize=fontsize)
#            ax3.imshow(w,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
            ax3.imshow(h,cmap=plt.cm.Blues)
            ax3.set_title(r'$h$', fontsize=fontsize)
#            im4=ax4.imshow(h,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
#            cbar_ax2 = fig.add_axes([0.85, 0.54, 0.03, 0.35])
#            fig.colorbar(im1, cax=cbar_ax2)
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()