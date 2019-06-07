# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:45:37 2016

For handling AUTO continuation for BWH model with tradeoff parameter chi

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
#=================
# AUTO tlm model
#=================
from auto import run,merge,rl,plot
import numpy as np
import deepdish.io as dd


def scan_for_betas(betas,name,resolution=11):
	data = []
	for beta in betas:
		chi_range,p_range = scan_chi_range(beta,resolution)
		data.append(p_range)
	data.insert(0,chi_range)
	data = np.array(data)
	np.savetxt(name,data)

def scan_chi_range(beta,resolution):
    chi_range = np.linspace(0,1,resolution)
    p_range = np.zeros(resolution)
    for i,chi in enumerate(chi_range):
        p_range[i]=obtain_Turing_prec_to_chi(chi,beta)
    return chi_range,p_range


def obtain_Turing_prec_to_chi(chi,beta=1.0,a=0,omegaf=1.0,HBnum=1,
                             NMX=2000,NPR=1000,DSMAX=0.0151,DS=0.0002,pmax=10.0,
                             model='bwh_to_EQK',parset='bwh_set4'):
    """ Loads parameters from parset hdf5 dictionary file, save in parameters txt file and perform continuation """
    p = loadparset(parset)
    p['chi']=chi
    p['a']=a
    p['omegaf']=omegaf
    p['beta']=beta
    a = conv_dict_to_parset_array(p)
    save_parset_to_txt_file(a)
    st0 = run(e=model,c=model,NMX=NMX,NPR=NPR,DS=DS,DSMAX=DSMAX,IPS=1,ISP=2)
    ind = 1
    while (len(st0("HB"))<HBnum) and (ind<100):
        if ind%2:
            print 'at %d!' %ind
            st0 = run(e=model,c=model,NMX=NMX,NPR=NPR,DS=DS,DSMAX=DSMAX/(1+ind/20.),IPS=1,ISP=2)
        else:
            print 'neg'
            st0 = run(e=model,c=model,NMX=NMX,NPR=NPR,DS=-DS,DSMAX=DSMAX/(1+ind/20.),IPS=1,ISP=2)
        ind+=1
    if st0("HB"):
        p = st0("HB1")['p']
    else:
        p = None
    return p


def GetTuring(name,step=100,HBnum=1,ICP1=['P','dummy','L','k'],DSMAX1=0.0173,NMX1=3000,UZZ=[]):
    st0 = run(e=name,c=name,NMX=NMX1,NPR=1000,DS=0.000001,DSMAX=DSMAX1,IPS=1,ISP=2,UZR=UZZ,ICP=ICP1)
    ind = 1
    while (len(st0("HB"))<HBnum) and (ind<100):
        if ind%2:
            print 'at %d!' %ind
            st0 = run(e=name,c=name,NMX=NMX1,NPR=1000,DS=0.000133,DSMAX=DSMAX1/(1+ind/20.),IPS=1,ISP=2,UZR=UZZ,ICP=ICP1)
        else:
            print 'neg'
            st0 = run(e=name,c=name,NMX=NMX1,NPR=1000,DS=-0.000133,DSMAX=DSMAX1/(1+ind/20.),IPS=1,ISP=2,UZR=UZZ,ICP=ICP1)
        ind+=1
    if st0("HB"):
        st1 = run(st0('HB%d'%min(HBnum,len(st0("HB")))),ISP=1,IPS=2,ICP=ICP1,NMX=step,NPR=1000,DSMAX=DSMAX1,DS=0.00001)
    else:
        st1 = []
    return st0,st1

def plotChi(chi=[0.0,0.5,1.0],
            model='bwh_to_EQK',parset='bwh_set4',
            pmax=10.0,DS=0.0002):
    chi02,ps2 = scanBif(chi=chi[0],parset=parset,pmax=pmax,model=model,DS=DS)
    chi05,ps5 = scanBif(chi=chi[1],parset=parset,pmax=pmax,model=model,DS=DS)
    chi07,ps7 = scanBif(chi=chi[2],parset=parset,pmax=pmax,model=model,DS=DS)
    plot(chi02+chi05+chi07)
    return [chi02,chi05,chi07],[ps2,ps5,ps7]

def scanBif(chi=0.0,a=0,omegaf=1.0,beta=1.0,
            NMX=2000,NPR=1000,DSMAX=0.0113,DS=0.0002,pmax=10.0,
            model='bwh_to_EQK',parset='bwh_set4'):
    """ Loads parameters from parset hdf5 dictionary file, save in parameters txt file and perform continuation """
    p = loadparset(parset)
    p['chi']=chi
    p['a']=a
    p['omegaf']=omegaf
    p['beta']=beta
    a = conv_dict_to_parset_array(p)
    save_parset_to_txt_file(a)
    bif = run(e=model,c=model,NMX=NMX,NPR=NPR,DSMAX=DSMAX,DS=DS,RL1 = pmax)
    return bif,p

def integrate(ic,p,chi=0.0,a=0,omegaf=1.0,NMX=10000,NPR=100,DSMAX=0.0113,DS=0.0002,Tmax=10000.0,model='tlm',parset='tlm_set10'):
    p = loadparset(parset)
    p['chi']=chi
    p['a']=a
    p['omegaf']=omegaf
    a = conv_dict_to_parset_array(p)
    save_parset_to_txt_file(a)
    result = run(e=model,c=model,NMX=NMX,NPR=NPR,DSMAX=DSMAX,DS=DS,RL1 = Tmax,ICP=[14],IPS=-2,U={'b':ic[0],'s1':ic[1],'s2':ic[2]})
    return result,ps




def saveBif(bif,fname="bwh_tf_bif",Ps=None):
    import deepdish as dd
    continuation=[]
    for i,branch in enumerate(bif):
        print "Branch",i+1, " of shape ", branch.toarray().shape
        LABs = []
        TYs = []
        for solution in branch:
            LABs.append(solution['LAB'])
            TYs.append(solution['TY name'])
        data = branch.toarray()
        branch_dict={'LAB':LABs,'TY':TYs,'p':data[0],'L2norm':data[1],
                      'b':data[2],'w':data[3],'h':data[4],
                      'sol':False}
        if Ps is not None:
            branch_dict['Ps']=Ps
        continuation.append(branch_dict)
    dd.io.save(fname+'.hdf5', continuation)

def saveHopf(hb,fname="tlm_hopf"):
    if type(hb)==dict:
        savePs=True
        Ps = hb['Ps']
        branch = hb['tur']
    else:
        savePs=False
        branch=hb
    import deepdish as dd
    turing = branch()
    p=[]
    t=[]
    b1sol=[]
    b2sol=[]
    s1sol=[]
    s2sol=[]
    b1=[]
    b2=[]
    s1=[]
    s2=[]
    nx = []
    lx = []
    for solution in turing:
        p.append(solution['p'])
        t.append(solution['t'])
        nx.append(len(solution['t']))
        lx.append(solution['L'])
        b1sol.append(solution['b1'])
        b2sol.append(solution['b2'])
        s1sol.append(solution['s1'])
        s2sol.append(solution['s2'])
        b1.append(np.mean(solution['b1']))
        b2.append(np.mean(solution['b2']))
        s1.append(np.mean(solution['s1']))
        s2.append(np.mean(solution['s2']))
    p = np.array(p)
    t  = np.array(t)
    nx = np.array(nx)
    lx = np.array(lx)
    b1sol=np.array(b1sol)
    b2sol=np.array(b2sol)
    s1sol=np.array(s1sol)
    s2sol=np.array(s2sol)
    b1=np.array(b1)
    b2=np.array(b2)
    s1=np.array(s1)
    s2=np.array(s2)
    solution = {'sol':True,'p':p,'t':t,
                'b1sol':b1sol,'b2sol':b2sol,'s1sol':s1sol,'s2sol':s2sol,
                'b1':b1,'b2':b2,'s1':s1,'s2':s2,
                'lx':lx,'nx':nx,'nd':1}
    if savePs:
        solution['Ps']=Ps
    dd.io.save(fname+'.hdf5', [solution])

def loadparset(parset):
    return dd.load(parset+'.hdf5')

def conv_dict_to_parset_array(ps):
    a = np.array([ps['p'],ps['lamb_max'],ps['lamb_min'],ps['eta'],
                  ps['nuw'],ps['nuh'],ps['rhow'],ps['rhoh'],
                  ps['gamma'],ps['alpha'],ps['f'],ps['q'],
                  ps['w_wp'],ps['w_fos'],ps['w_fc'],
                  ps['mu_s_max'],ps['omegaf'],
                  ps['chi'],ps['beta'],ps['a'],
                  ps['dw'],ps['dh'],ps['del_to']])
    return a
def save_parset_to_txt_file(a):
    np.savetxt("bwh_tf_parameters.txt",a.T)
