#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 17:03:23 2018

@author: ohm
"""

def setEs(args):
    Es={'rhs':args.rhs,
        'n':tuple(args.n),
        'l':tuple(args.l),
        'bc':args.bc,
        'it':args.it,
        'dt':args.dt,
        'analyze':True,
        'verbose':args.verbose,
        'setPDE':True}
    return Es
    

def main(args):
    from bwhModel import bwhModel
    Es=setEs(args)
    Ps=args.Ps
    Vs=args.Vs
    m = bwhModel(Ps,Es,Vs)
    step = args.step*m.p['conv_T_to_t']
    max_time = args.finish*m.p['conv_T_to_t']
    savefile=args.fname
    check_convergence=args.check_convergence
    create_movie=args.movie
    sol = m.integrate(step=step,max_time=max_time,savefile=savefile,
                      create_movie=create_movie,
                      check_convergence=check_convergence,
                      p=args.p,chi=args.chi)
    if args.fname is not None:
        m.save_state(sol,args.fname)
    return 0

def add_parser_arguments(parser):
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Turn on debuging messages")
    parser.add_argument("--check_convergence",
                        action="store_true",
                        dest="check_convergence",
                        default=False,
                        help="Run pde time integration")
    parser.add_argument("-f", "--fname",
                        type=str, nargs='?',
                        dest="fname",
                        default=None,
                        help="Save bwh integration in fname")
    parser.add_argument("-n", 
                        nargs='+', type=tuple,
                        dest="n",
                        default=(256,),
                        help="Number of grid points")
    parser.add_argument("-l",
                        nargs='+', type=tuple,
                        dest="l",
                        default=(128.0,),
                        help="Systm length")
    parser.add_argument("--Vs",
                        type=str, nargs='?',
                        dest="Vs",
                        default="random",
                        help="Initial condition name of integration")
    parser.add_argument("--Ps",
                        type=str, nargs='?',
                        dest="Ps",
                        default="auto/bwh_set2.hdf5",
                        help="Initial condition name of integration")
    parser.add_argument("--rhs",
                        type=str, nargs='?',
                        dest="rhs",
                        default="oz_EQK",
                        help="RHS type of the model")
    parser.add_argument("--bc",
                        type=str, nargs='?',
                        dest="bc",
                        default="periodic",
                        help="Boundary conditions")
    parser.add_argument("--it",
                        type=str, nargs='?',
                        dest="it",
                        default="pseudo_spectral",
                        help="Integrator type")
    parser.add_argument("--movie",
                        action="store_true",
                        dest="movie",
                        default=False,
                        help="Make a movie")
    parser.add_argument('--dt',
                        dest='dt',
                        type=float,
                        default=0.1,
                        help='Time step as multiplication of Von Neumann')
    parser.add_argument('-p','--prec',
                        dest='p',
                        type=float,
                        default=0.1,
                        help='Precipitation')
    parser.add_argument('-c','--chi',
                        dest='chi',
                        type=float,
                        default=0.0,
                        help='Trade-off parameter')
    parser.add_argument('--finish',
                        dest='finish',
                        type=float,
                        default=100.0,
                        help='End time in years')
    parser.add_argument('--step',
                        dest='step',
                        type=float,
                        default=0.1,
                        help='Ouput solution at given interval in years')
    parser.add_argument('--start',
                        dest='start',
                        type=float,
                        default=0.0,
                        help='Ouput solution at given interval')

    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    main(parser.parse_args())