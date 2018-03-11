# set3 - No Eta
# Set4 - With Eta=7.0
# Set6 - Reasonable values for minimum P value for growth, no imaginary eigenvalues
# Shallow Z1 and Z2

#from math import sqrt
Z1 = 5.0
Z2 = 400.0
npor = 0.55
s_wp = 0.065 # wilting point
s_fos = 0.175 # fully open-stomata point
s_h  = 0.05 # residual water content / hygroscopic point
M0    = 20.0
mu_s_max = 0.30
Theta_fc = 0.16
s_fc = Theta_fc/npor
Lambda = M0/(s_wp*npor*Z2)

dimpar = {'Lambda_min' : 0.35,
          'Lambda_max' : 0.40,
          'K'      : 0.40,
          'E'      : 5.0,
          'M0'      : M0,
          'N'      : 100.0,
          'DB'     : 1.2,
          'DW'     : 150.0,
          'DH'     : 150.0,
          'R'      : 0.5,
          'Gamma'  : 4.0,
          'P'      : 10.0,
          'A'      : 10.0,
          'Q'      : 10.0,
          'F'      : 10.0,
          'nZ1'    : npor*Z1,
          'nZ2'    : npor*Z2,
          'Ks'     : 5000.0,
          's_wp'   : s_wp,
          'seasonality' : 0.0,
          's_fos'  : s_fos,
          's_fc'   : s_fc,
          'mu_s_max' : mu_s_max,
          's_h'    : s_h,
         }


#print dimpar
# Calculation of non dimensional parameters from the dimensional parameters
def update_par():
    par = {'lamb_max' : (dimpar['Lambda_max']*dimpar['nZ2'])/dimpar['M0'],
           'lamb_min' : (dimpar['Lambda_min']*dimpar['nZ2'])/dimpar['M0'],
           'eta'      : dimpar['E']*dimpar['K'],
           'p'        : dimpar['P']/(dimpar['nZ1']*dimpar['M0']),
           'nu'       : dimpar['N']/dimpar['M0'],
           'rho'      : dimpar['R'],
           'kappa'    : dimpar['Ks']/(dimpar['nZ1']*dimpar['M0']),
           'c'        : 2,
           'alpha'    : 10.0,
           'q'        : 10.0,
           'f'        : 10.0,
           'zeta'     : Z1/Z2,
           'gamma'    : (dimpar['Gamma']*dimpar['K'])/dimpar['M0'],
           'dw'       : dimpar['DW']/dimpar['DB'],
           'dh'       : dimpar['DH']/dimpar['DB'],
           'conv_P'   : (dimpar['nZ1']*dimpar['M0']),
           'conv_T_to_t'   : dimpar['M0'],
           's_wp'     : dimpar['s_wp'],
           's_fos'    : dimpar['s_fos'],
           's_fc'     : dimpar['s_fc'],
           's_h'      : dimpar['s_h'],
           'mu_s_max' : mu_s_max,
           'chi'      : 0.0,
           'beta'     : 1.0,
           'omegaf'   : 1.0,
           'a'        : dimpar['seasonality'],
           'dimpar'   :dimpar
           }
    return par

par=update_par()


def savecsvdimpar(fname,dimpar):
    import csv
    with open(fname+'.csv', 'wb') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['Parameter','value'])
        writer.writeheader()
        for k in dimpar:
            writer.writerow({'Parameter': k , 'value': dimpar[k]})
    
def saveParmSet(fname,par,text=None,saveroot=False):
    from deepdish.io import save
    if saveroot:
        from tlmModel import tlmModel,Es_normal
        m = tlmModel(Es=Es_normal,Ps=par,Vs=None)
        t,result = m.ode_integrate([0.9,0.2,0.2])
        b1,s1,s2 = result.T[-1]
        par['b1']=b1
        par['s1']=s1
        par['s2']=s2
    if text is not None:
        par['text']=text
    save("./auto/"+fname+".hdf5",par)
def loadParmSet(fname):
    from deepdish.io import load
    return load(fname)

if __name__ == '__main__':
    par=update_par()
    print "Nondimensional:",
    print par
    saveParmSet('bwh_set1',par,saveroot=False)
#    import numpy as np
#    p=par
#    a = np.array([p['lamb_max'],p['lamb_min'],p['eta'],
#                  p['p'],p['nu'],p['rho'],p['kappa'],
#                  p['c'],p['zeta'],p['gamma'],
#                  p['s_wp'],p['s_fos'],p['s_fc'],
#                  p['s_h'],p['mu_s_max'],
#                  p['chi'],p['beta'],p['iota'],
#                  p['omegaf'],p['delta_s']])
#    np.savetxt("./auto/tlm_parameters.txt",a.T)
##    savemat("./p2p/b2s2_parameters.mat",p)