#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 07:08:40 2018

@author: ohm
"""
from bwhModel import bwhModel
import numpy as np
set2 = 'auto/bwh_set2.hdf5'
def simdrought(prec_i,prec_f,delta_p,delta_year,chi,
               Ps=set2,
               n=(512,512),l=(256.0,256.0),
               Vs_initial="random",rhs="oz_EQK",
               bc="periodic",it="pseudo_spectral",
               first_time = 100.0,tol=1.0e-8,add_noise=0.01,
               fname="cont",verbose=True,
               savefile=None,create_movie=False,
               send_email=None):
    import deepdish.io as dd
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    Es={'rhs':rhs,'n':n,'l':l,'bc':bc,'it':it,
        'dt':0.1,'verbose':verbose,'analyze':False,'setPDE':True}
    if type(Vs_initial)==str:
        fname = fname+"_"+Vs_initial
    prec_gradient_down = np.arange(prec_i,prec_f-delta_p,-delta_p)
    time_span = np.arange(delta_year,len(prec_gradient_down)*delta_year+delta_year,delta_year)
    m = bwhModel(Vs=Vs_initial,Es=Es,Ps=Ps)
    m.setup['verbose']=verbose
    yr=m.p['conv_T_to_t']
    # Converging on the first solution using integration and then root
    Vs_init = m.integrate(m.initial_state,check_convergence=True,
                          max_time=first_time*yr,p=prec_i,chi=chi)
#    Es['rhs']="oz_EQK_relax"
#    m = bwhModel(Vs=Vs_initial,Es=Es,Ps=Ps)
#    m.setup['verbose']=verbose
    Vs = Vs_init.copy()
    b_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    w_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    h_sol = np.zeros((len(prec_gradient_down),n[0],n[1]))
    if create_movie and savefile is not None:
        savefile_base=savefile
    for i,prec in enumerate(prec_gradient_down):
        print "Integration for p =",prec
        if create_movie and savefile is not None:
            savefile=savefile_base+"_p{:4.3f}".format(prec).replace(".","_")
        b,w,h=m.split_state(Vs)
        if add_noise is not None:
            b=b+add_noise*np.random.random(size=b.shape)
            w=w+add_noise*np.random.random(size=w.shape)
#            h=h+add_noise*np.random.random(size=h.shape)
        Vs = np.ravel((b,w,h))
        Vs_new=m.integrate(initial_state=Vs,
                           max_time=delta_year*yr,step=yr,
                           check_convergence=False,
                           savefile=savefile,create_movie=False,
                           p=prec,chi=chi)
#        if m.converged_relaxation==False:
#            time,result=m.pseudo_spectral_integrate(initial_state=Vs,
#                                                    finish=delta_year*yr,
#                                                    step=yr,
#                                                    p=prec,chi=chi)
#            Vs_new=result[-1]
        b,w,h=m.split_state(Vs_new)
        b_sol[i]=b
        w_sol[i]=w
        h_sol[i]=h
    dd.save(fname+".hdf5",{'p':prec_gradient_down,"T":time_span,'chi':chi,
                           'Ps_dimensional':m.p['dimpar'],
                           'n':n,'l':l,
                           'b':b_sol,
                           'w':w_sol,
                           'h':h_sol},
            compression='blosc')
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Drought bwh simulations finished",
                    "",
                    "Drought bwh simulations finished and saved to:", fname
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def main(args):
    simdrought(args.prec_i,args.prec_f,args.delta_p,
               args.delta_year,args.chi,add_noise=args.noise,
               create_movie=args.create_movie,savefile=args.savefile,
               fname=args.fname,verbose=args.verbose,send_email=args.send_email)
    return 0

def add_parser_arguments(parser):
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Turn on debuging messages")
    parser.add_argument("--create_movie",
                        action="store_true",
                        dest="create_movie",
                        default=False,
                        help="Create movie at each simulation")
#    parser.add_argument("--check_convergence",
#                        action="store_true",
#                        dest="check_convergence",
#                        default=False,
#                        help="Run pde time integration")
    parser.add_argument("-f", "--fname",
                        type=str, nargs='?',
                        dest="fname",
                        default="bwh_sim_drought",
                        help="Save bwh integration in fname")
#    parser.add_argument("-n", 
#                        nargs='+', type=tuple,
#                        dest="n",
#                        default=(256,),
#                        help="Number of grid points")
#    parser.add_argument("-l",
#                        nargs='+', type=tuple,
#                        dest="l",
#                        default=(128.0,),
#                        help="Systm length")
#    parser.add_argument("--Vs",
#                        type=str, nargs='?',
#                        dest="Vs",
#                        default="random",
#                        help="Initial condition name of integration")
#    parser.add_argument("--Ps",
#                        type=str, nargs='?',
#                        dest="Ps",
#                        default="auto/bwh_set2.hdf5",
#                        help="Initial condition name of integration")
#    parser.add_argument("--rhs",
#                        type=str, nargs='?',
#                        dest="rhs",
#                        default="oz_EQK",
#                        help="RHS type of the model")
#    parser.add_argument("--bc",
#                        type=str, nargs='?',
#                        dest="bc",
#                        default="periodic",
#                        help="Boundary conditions")
#    parser.add_argument("--it",
#                        type=str, nargs='?',
#                        dest="it",
#                        default="pseudo_spectral",
#                        help="Integrator type")
#    parser.add_argument("--movie",
#                        action="store_true",
#                        dest="movie",
#                        default=False,
#                        help="Make a movie")
    #prec_i,prec_f,years,delta_year,chi
    parser.add_argument("--send_email",
                        type=str, nargs='?',
                        dest="send_email",
                        default=None,
                        help="Password to send email notice")
    parser.add_argument("--savefile",
                        type=str, nargs='?',
                        dest="savefile",
                        default=None,
                        help="Save each simulation to base name")
    parser.add_argument('--chi',
                        dest='chi',
                        type=float,
                        default=0.0,
                        help='Trade-off parameter value')
    parser.add_argument('--noise',
                        dest='noise',
                        type=float,
                        default=0.01,
                        help='Noise addition in every step')
    parser.add_argument('--prec_i',
                        dest='prec_i',
                        type=float,
                        default=3.5,
                        help='Initial value for precipitation')
    parser.add_argument('--prec_f',
                        dest='prec_f',
                        type=float,
                        default=1.0,
                        help='Final value for precipitation')
    parser.add_argument('--delta_p',
                        dest='delta_p',
                        type=float,
                        default=0.125,
                        help='Number of years')
    parser.add_argument('--delta_year',
                        dest='delta_year',
                        type=int,
                        default=40,
                        help='Number of years')
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    main(parser.parse_args())