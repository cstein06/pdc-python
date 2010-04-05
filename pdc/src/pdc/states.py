# -*- coding:utf-8 -*-
"""
Created on Apr 5, 2010

@author: Stein
"""

from numpy import *

import time

import matplotlib.pyplot as pp

import scipy.stats as st
from scipy.linalg.basic import det
from scipy.io import savemat

import pdc.plotting as pl
import pdc.ar_fit as fit_
from pdc.globals import *
from pdc.ar_data import ar_data 
import pdc.analysis as an_


def window_analysis(data, **args):
    
    pr_.do_log = False 
    
    read_args(args)
    
    print '\nTotal data shape:', data.shape
    n, T = data.shape
    
    win = int(pr_.window_size*pr_.sample_f)
    
    nwins = T/win
    
    print '\nCalculating first window verbosely (others will be silent):\n'
    aux_v = pr_.v
    
    aux = an_.measure(data[:,0*win:0*win+win], ss = False)
    
    resp = []
    if type(aux) is 'tuple':
        for i in arange(len(aux)):
            resp.append(zeros([nwins]+list(aux[i].shape), dtype = aux[i].dtype))
            resp[i][0] = aux
    else:
        resp = [zeros([nwins]+list(aux.shape), dtype = aux.dtype)]
            
    for i in arange(0,nwins):
        
        #print i*win
        
        aux = an_.measure(data[:,i*win:i*win+win], ss = False)
        for j in arange(len(resp)):
            resp[j][i] = aux[j]
        
        print '\nProcessed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
        
        pr_.v = False
    
    pr_.v = aux_v
    
    return resp

def mean_states(mes, states):
        
    nwins,n,n,nf = mes.shape
    
    if mes.dtype == 'complex':
        print 'Taking mean of complex measure!`'
    
    mpdc = zeros([len(pr_.valid_states), n, n, nf], dtype = mes.dtype)
    spdc = zeros([len(pr_.valid_states), n, n, nf], dtype = mes.dtype)
    
    states = states[:nwins]
    
    if (nwins > size(states)):
        print 'More windows than states in the states file!'
    
    for i in pr_.valid_states:
        if sum(states == i) > 0:
            auxi = pr_.st_dict[i]
            mpdc[auxi] = mean(mes[states == i], 0)
            spdc[auxi] = std(mes[states == i], 0)
            
    return mpdc, spdc

def mean_states_list(mes, estag, maxe = 6):
        
    mpdc = []
    spdc = []
    for i in arange(len(mes)):
        aux = mean_states(mes[i], estag)
        mpdc.append(aux[0])
        spdc.append(aux[1])
            
    return mpdc, spdc
        
def states_analysis(data, states, **args):
        
    read_args(args)

    tim = time.clock()
    
    
    result = window_analysis(data)
    
    nwins = result[0].shape[0]
    
    states = states[:nwins]
    
    if pr_.valid_states is None:
        set_params(valid_states = list(sorted(set(states[states >= 0]))))
    
    nstates = zeros(len(pr_.valid_states))
    
    pr_.st_dict = {}
    for i in arange(len(pr_.valid_states)):
        pr_.st_dict[pr_.valid_states[i]] = i
        nstates[i] = sum(states == pr_.valid_states[i])
    
    #nstates = histogram(states, bins=arange(1,8))[0]
    
    print '\nNumber of windows used:', states.shape[0]
    
    mpdc, spdc = mean_states_list(result, states)
    
    #nf = mpdc.shape[3]
    
    print '\nNumber of windows per state:', nstates
        
    if pr_.do_plot:
        if pr_.plot_states is None:
            pr_.plot_states = pr_.valid_states
        for st in pr_.plot_states:
            i = pr_.st_dict[st]
            if nstates[i] > 0:
                pr_.plot_color = pr_.state_colors[i]
                res_.mes = mpdc[0][i]
                pl.plot_all()
        pp.show()
        
    if len(mpdc) == 1:
        result = result[0]
        mpdc = mpdc[0]
        spdc = spdc[0]
                
    print '\nTotal time in secs:', time.clock() - tim
    
    if pr_.do_states_log:
        log_windows_results(result, mpdc, spdc, nstates)
    
    return result, mpdc, spdc, nstates
