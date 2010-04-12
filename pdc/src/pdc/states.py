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
    
    aux = an_.measure(data[:,0*win:0*win+win])
    if type(aux) is not tuple:
        aux = [aux]
    
    resp = []
    for i in arange(len(aux)):
        resp.append(zeros([nwins]+list(aux[i].shape), dtype = aux[i].dtype))
        resp[i][0] = aux[i]
            
    for i in arange(0,nwins):
        
        pr_.v = False
        #print i*win
        
        aux = an_.measure(data[:,i*win:i*win+win])
        if type(aux) is not tuple:
            aux = [aux]
            
                  
        for j in arange(len(resp)):
            resp[j][i] = aux[j]
            
        if pr_.ss is True:
            for k in arange(n):
                resp[0][i,k,k,:] = res_.ss[k,k,:]
        
        print '\nProcessed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
        
    
    pr_.v = aux_v
    
    return resp

def mean_states(mes, states):
    
    nwins = mes.shape[0]
        
    if mes.dtype == 'complex':
        print 'Taking mean of complex measure!`'
    
    mpdc = zeros([len(pr_.valid_states)]+list(mes[0].shape), dtype = mes.dtype)
    spdc = zeros([len(pr_.valid_states)]+list(mes[0].shape), dtype = mes.dtype)
    
    states = states[:nwins]
    
    if (nwins > size(states)):
        print 'More windows than states in the states file!'
    
    for i in arange(len(pr_.valid_states)):
        if sum(states == i) > 0:
            #auxi = pr_.st_dict[pr_.valid_states[i]]
            auxi = i
            mpdc[auxi] = mean(mes[states == pr_.valid_states[i]], 0)
            spdc[auxi] = std(mes[states == pr_.valid_states[i]], 0)
            
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
                res_.ss = res_.mes
                pr_.plot_color = pr_.state_colors[i]
                res_.mes = mpdc[0][i]
                pl.plot_all()
        pp.show()
        
#    if len(mpdc) == 1:
#        result = result[0]
#        mpdc = mpdc[0]
#        spdc = spdc[0]
                
    print '\nTotal time in secs:', time.clock() - tim
    
    if pr_.do_states_log:
        log_windows_results(result, mpdc, spdc, nstates, states)
    
    return result, mpdc, spdc, nstates


def states_analysis_bind(data, states, **args):
        
    read_args(args)

    tim = time.clock()
    
    n, nd = data.shape
    win = int(pr_.window_size*pr_.sample_f)
    
    states = states[:int(nd/win)]
    
    if pr_.valid_states is None:
        set_params(valid_states = list(sorted(set(states[states >= 0]))))
    
    nstates = zeros(len(pr_.valid_states))
    pr_.st_dict = {}
    for i in arange(len(pr_.valid_states)):
        pr_.st_dict[pr_.valid_states[i]] = i
        nstates[i] = sum(states == pr_.valid_states[i])
        
    
    mpdc = zeros([len(pr_.valid_states), n, n, pr_.nf]) 
    spdc = zeros([len(pr_.valid_states), n, n, pr_.nf]) 
    for i in arange(len(pr_.valid_states)):
        if nstates[i] > 0:  
            datast = array([]).reshape(n,-1)
            for j in where(states == pr_.valid_states[i])[0]:
                datast = concatenate((datast, data[:,win*j:win*j+win]), axis = 1)
            
            mpdc[i], th, ic1, ic2 = an_.measure_full(datast)
            spdc[i] = mpdc[i] - ic1 
    
    result = []
    
    print '\nNumber of windows per state:', nstates
        
#    if pr_.do_plot:
#        if pr_.plot_states is None:
#            pr_.plot_states = pr_.valid_states
#        for st in pr_.plot_states:
#            i = pr_.st_dict[st]
#            if nstates[i] > 0:
#                pr_.plot_color = pr_.state_colors[i]
#                res_.mes = mpdc[0][i]
#                pl.plot_all()
    pp.show()
                
    print '\nTotal time in secs:', time.clock() - tim
    
    if pr_.do_states_log:
        log_windows_results([result], [mpdc], [spdc], nstates, bind = True)
    
    return result, mpdc, spdc, nstates
