# -*- coding:utf-8 -*-
"""
Created on 06/11/2009

@author: Carlos Stein
"""

__all__ = ['event_adaptative', 'adaptative_ar']

from numpy import *
import scipy.signal as sig
import time
from matplotlib.pyplot import gcf
import matplotlib.pyplot as pp

from pdc import *

def event_adaptative(data, events, se = 300, preproc = False, 
                     win = (-10,200), **args):
    read_args(args)
    
    n,nd = data.shape
    
    ne = events.size
    
    nw = win[1]-win[0]+1
    
    datan = zeros([ne, n, nw])
    for i in arange(ne):
        datan[i] = data[:,events[i]+win[0]:events[i]+win[1]+1]
    
    return adaptative_ar(datan, se = se, preproc = preproc)
    
def adaptative_ar(data, se = 100, preproc = False, **args):
    '''data(m,n,nd) -> data, m = #trials, n = #channels, nd = #time samples
       se -> efective sample memory of adaptative model'''
    
    read_args(args)
    
    m,n,nd = data.shape
    
    if preproc:
        data = pre_proc_ding_99(data)
    
    A, er = adapt_ar(data, pr_.maxp, se)
    
#     Plot some samples
#    for i in arange(pr_.maxp, nd, step):
#        res_.mes = vars()[pr_.alg + '_alg'](A[i], er[i], metric = 'diag')
#        pr_.power = False
#        plot_all()
#        canvas = gcf().canvas
#        canvas.start_event_loop(timeout=0.2)
#
#        pr_.v = False
        #time.sleep(1)
    
#    pr_.v = True
        
    resg = zeros([nd,n,n,pr_.nf], dtype = 'complex')
    
    for i in arange(pr_.maxp, nd):
        pr_.power = False
        resg[i] = globals()[pr_.alg + '_alg'](A[i], er[i], metric = pr_.metric)
    
    pp.figure()
    plot_coherogram(resg)
    pp.show()
    
    return A, er
         

def pre_proc_ding_99(data):
    
    m,n,nd = data.shape
    
    #normalize per trial
    data = sig.detrend(data)
    data = data/std(data, axis = 2).reshape(m,n,1)
    
    #normalize per time
    if m > 1:
        data = data - mean(data, 0).reshape(1,n,nd)
        data = data/std(data, axis = 0).reshape(1,n,nd)
    
    return data
    
    
            
