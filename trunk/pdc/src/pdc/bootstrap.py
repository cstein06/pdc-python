# -*- coding:utf-8 -*-
"""
Created on 22/10/2009

@author: Carlos Stein
"""

from numpy import *
import time
from pdc.ar_fit import nstrand

from pdc.ar_data import ar_data

def bootstrap(method_func, nd, nm, A, er, 
              nf, alpha = 0.05, metric = None):
    '''Estatistica bootstrap'''
    
    n = A.shape[0]
    mes = empty([nm, n, n, nf])
    #th  = empty([nm, n, n, nf])
    #ic1 = empty([nm, n, n, nf])
    #ic2 = empty([nm, n, n, nf])
    maxp = A.shape[2]
    tbegin = time.clock()
    for i in range(nm):
        if (i%100 == 0):
            print 'nm:', i, 'time:', time.clock()-tbegin
        #Generate data from AR
        data = ar_data(A, er, nd)
        #Estimate AR parameters with Nuttall-Strand
        Aest, erest = nstrand(data, maxp = maxp)
        #Calculate the connectivity and statistics
        if (metric == None):
            mes[i] = abs(method_func(Aest, erest, nf = nf))**2
        else:
            mes[i] = abs(method_func(Aest, erest, nf = nf, metric = metric))**2

    so = sort(mes, axis = 0)
    ic1 = so[(alpha/2)*nm]
    ic2 = so[(1-alpha/2)*nm]
    #bvar = var(mes, axis = 0)

    if (metric == None):
        mes0 = abs(method_func(A, er, nf = nf))**2
    else:
        mes0 = abs(method_func(A, er, nf = nf, metric = metric))**2
    
    th = ones(ic1.shape)*0#0.5/sqrt(nd) #TODO: pensar melhor como fazer.
    
    return mes0, th, ic1, ic2
