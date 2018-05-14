# set1 - from OZ paper parameters

dimpar = {'Lambda_min' : 0.025,
          'Lambda_max' : 0.030,
          'K'      : 2.0/3.0,
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
          'P'      : 100.0,
          'A'      : 120.0,
          'Q'      : 1.2,
          'F'      : 0.01,
          's_wp'   : 0.085,
          's_fos'  : 0.175,
          'mu_s_max' : 0.30,
          'seasonality': 0.0,
          'npor'   : 0.30,
          'del_to' : 0.3,
          'Z'      : 1000.0,
         }
dimpar['s_fc'] = 0.16/dimpar['npor'] # Theta_fc/npor
dimpar['nZ']=  dimpar['npor']*dimpar['Z']
#print dimpar
# Calculation of non dimensional parameters from the dimensional parameters
def update_par():
    par = {'lamb_max' : 1.0,
           'lamb_min' : dimpar['Lambda_min']/dimpar['Lambda_max'],
           'eta'      : dimpar['E']*dimpar['K'],
           'del_to' : dimpar['del_to'],
           'p'        : (dimpar['P']*dimpar['Lambda_max'])/(dimpar['M0']**2),
           'nuw'      : dimpar['NW']/dimpar['M0'],
           'nuh'      : dimpar['NH']/dimpar['M0'],
           'rhow'     : dimpar['RW'],
           'rhoh'     : dimpar['RH'],
           'alpha'    : dimpar['A'] / dimpar['M0'],
           'q'        : dimpar['Q'] / dimpar['K'],
           'f'        : dimpar['F'],
           'gamma'    : (dimpar['Gamma']*dimpar['K'])/dimpar['M0'],
           'dw'       : dimpar['DW']/dimpar['DB'],
           'conv_T_to_t'   : dimpar['M0'],
           'dh'       : (dimpar['DH'] * dimpar['M0']) / (dimpar['DB'] * dimpar['Lambda_max']),
           'w_wp'     : dimpar['s_wp']*(dimpar['nZ']*dimpar['Lambda_max'])/dimpar['M0'],
           'w_fos'    : dimpar['s_fos']*(dimpar['nZ']*dimpar['Lambda_max'])/dimpar['M0'],
           'w_fc'     : dimpar['s_fc']*(dimpar['nZ']*dimpar['Lambda_max'])/dimpar['M0'],
           'mu_s_max' : dimpar['mu_s_max'],
           'chi'      : 0.5,
           'beta'     : 1.0,
           'omegaf'   : 1.0,
           'a'        : dimpar['seasonality'],
           'dimpar'   : dimpar
           }
    par['pc']=par['nuw']*((par['alpha']*par['f'])+par['nuh'])/(par['alpha']*par['f'])
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
    from scipy.io import savemat
    par=update_par()
    name = 'bwh_set2'
    par['pc']=par['nuw'] * ( (par['alpha'] * par['f']) + par['nuh'])/ (par['alpha'] * par['f'])
    print "Nondimensional:",
    print par
    saveParmSet(name,par,saveroot=False)
    savemat("./p2p/"+name+".mat",par)