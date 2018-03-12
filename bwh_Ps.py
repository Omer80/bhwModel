# set1 - from OZ paper parameters

dimpar = {'Lambda_min' : 0.025,
          'Lambda_max' : 0.030,
          'K'      : 0.666,
          'E'      : 1.5,
          'M0'     : 2.0,
          'NW'     : 1.5,
          'NH'     : 4.5,
          'DB'     : 0.1,
          'DW'     : 2.5,
          'DH'     : 4.0,
          'RW'     : 0.3,
          'RH'     : 0.8,
          'Gamma'  : 14.0,
          'P'      : 10.0,
          'A'      : 120.0,
          'Q'      : 1.2,
          'F'      : 0.01,
          'Ks'     : 5000.0,
          's_wp'   : 0.065,
          's_fos'  : 0.175,
          's_fc'   : 0.16/0.55, # Theta_fc/npor
          'mu_s_max' : 0.30,
          'seasonality': 0.0,
          'npor'   : 0.55,
          'Z'      : 1000.0,
         }

dimpar['nZ']=  dimpar['npor']*dimpar['Z']
#print dimpar
# Calculation of non dimensional parameters from the dimensional parameters
def update_par():
    par = {'lamb_max' : (dimpar['Lambda_max']*dimpar['nZ'])/dimpar['M0'],
           'lamb_min' : (dimpar['Lambda_min']*dimpar['nZ'])/dimpar['M0'],
           'eta'      : dimpar['E']*dimpar['K'],
           'p'        : dimpar['P']/(dimpar['nZ']*dimpar['M0']),
           'nuw'      : dimpar['NW']/dimpar['M0'],
           'nuh'      : dimpar['NH']/dimpar['M0'],
           'rhow'     : dimpar['RW'],
           'rhoh'     : dimpar['RH'],
           'kappa'    : dimpar['Ks']/(dimpar['nZ']*dimpar['M0']),
           'c'        : 2,
           'alpha'    : dimpar['A'] / dimpar['M0'],
           'q'        : dimpar['Q'] / dimpar['K'],
           'f'        : dimpar['F'],
           'gamma'    : (dimpar['Gamma']*dimpar['K'])/dimpar['M0'],
           'dw'       : dimpar['DW']/dimpar['DB'],
           'dh'       : (dimpar['DH'] * dimpar['M0']) / (dimpar['DB'] * dimpar['Lambda_max']),
           'conv_P'   : (dimpar['nZ']*dimpar['M0']),
           'conv_T_to_t'   : dimpar['M0'],
           's_wp'     : dimpar['s_wp'],
           's_fos'    : dimpar['s_fos'],
           's_fc'     : dimpar['s_fc'],
           'mu_s_max' : dimpar['mu_s_max'],
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
        b,w,h = result.T[-1]
        par['b']=b
        par['w']=w
        par['h']=h
    if text is not None:
        par['text']=text
    save("./auto/"+fname+".hdf5",par)
def loadParmSet(fname):
    from deepdish.io import load
    return load(fname)

if __name__ == '__main__':
    par=update_par()
    par['pc']=par['nuw'] * ( (par['alpha'] * par['f']) + par['nuh'])/ (par['alpha'] * par['f'])
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