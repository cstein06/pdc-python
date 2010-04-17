# -*- coding:utf-8 -*-
"""
This is a file of usage reference for the user

Mini instructions:

Use pr_ structure to set parameters. 
See available parameters in params.py->Params and their defaults.
Call measure_full(data)

Please refer to the manual to further explanations on usage.
"""

from pdc import *

from numpy import *

import matplotlib.pyplot as pp

def simple_simulation_analysis():
    '''Simple test of connectivity routines
    
    Here we simulate some data using the ar_data module, from a given
    MVAR matrix A and covariance matrix er.
    
    '''
    
    #Definition of the MVAR model
    A = array([[[0.2, 0],[0.3,-0.2],[0.3,-0.2]], 
               [[0, 0],[0.8,-0.1],[0.4,-0.1]],
               [[0, 0],[0.3,0.2],[0.4,0.1]]], dtype = float) 
    er = identity(3)
    
    #number of data points generated
    nd = 2000
    
    #number of frequency points analyzed
    pr_.nf = 64
    
    #error of the confidence interval and threshold
    pr_.alpha = 0.05
    
    #model order parameters
    n = A.shape[0]
    pr_.maxp = A.shape[2]
    
    #type of PDC used (refer to manual to see what it means)
    pr_.metric = 'diag'
    
    # Other parameters...
    pr_.ss = True
    pr_.plot_ic = True
    pr_.do_plot = True
    pr_.stat = 'asymp'
    
    # Set which connectivity method to use
    pr_.alg = 'coh'
    
    #Generate data from AR
    data = ar_data(A, er, nd)
    
    #Call any connectivity routine.
    measure_full(data)
        
    #If running file, to keep plot open.
    pp.show()
    
def simple_with_file():
    '''Simple test of connectivity routines
    
    Load file and calculate...
    '''
    
    file = 'your_file_here'
    
    data = loadtxt(file).T # Please use channel in rows, time in colums.
    
    #if it is a .mat file:
    #from scipy.io import loadmat
    #data = loadmat(file)['name_of_the_variable'].T
    
    # Parameters...
    pr_.maxp = 25
    pr_.alg = 'pdc'
    
    #Call any connectivity routine.
    measure_full(data)
        
    pp.show()
    
    
def window_simulation_analysis():
    
    root_dir = 'your_root_dir_for_data'
    pr_.stinput = 'data_dile'
    
    pr_.alg = 'coh'
    
    pr_.window_size = 10 # in seconds
    pr_.sample_f = 500 # in Hz
    
    pr_.nf = 256
    
    pr_.plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    pr_.plotf = 150
    pr_.plota = True
    pr_.do_window_log = True
    
    pr_.ordem_max = 25
    pr_.fixp = True # Fix model order for a lot of speed performance.
                    # To find which order is reasonable, run without it in some data and check order.

    data = loadtxt(root_dir+pr_.stinput).T
    
    print 'Data loaded:', pr_.stinput
    
    res = window_analysis(data)
    
    plot_coherogram(res)
    
    pp.show()
    
    return res
    
def states_simulation_analysis():
    
    root_dir = 'your_root_dir_for_data'
        
    pr_.stinput = 'data_dile'
    
    instates = 'states_file'
    
    pr_.alg = 'coh'
    
    pr_.window_size = 10
    pr_.nf = 250
    pr_.sample_f = 500
    
    pr_.plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    pr_.plot_states = array([1,2,3,4,5,6])
    pr_.plotf = 150
    pr_.plota = True
    pr_.do_window_log = True
    
    pr_.valid_states = [1,2,3,4,5,6]
    pr_.ordem_max = 25
    pr_.fixp = True
    
    pr_.metric = 'diag'
    
    data = loadtxt(root_dir+pr_.stinput).T
    states = loadtxt(root_dir+instates)
    
    print 'Data loaded:', pr_.stinput
    
    res, mea, stds, nstates = states_analysis(data, states)
    
    pp.figure()
    
    plot_coherogram(res, states)
    
    pp.show()
    
    return res, mea, stds, nstates


if __name__ == '__main__':
    
    simple_simulation_analysis()
    #states_simulation_analysis()
    #window_simulation_analysis()
    
    
