# -*- coding:utf-8 -*-
"""
This is a file you user could use for your code.

Mini instructions:

Use an_.set_params() to set parameters. 
See available parameters in globals.py->Params.
Call an_.measure_full(data)

"""

import pdc.analysis as an_
import pdc.ar_data as ar_
import pdc.ar_fit as fit_
import pdc.plotting as pl_
import pdc.asymp as ass_

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
    nf = 64
    
    #error of the confidence interval and threshold
    alpha = 0.05
    
    #model order parameters
    n = A.shape[0]
    maxp = A.shape[2]
    
    #type of PDC used (refer to manual to see what it means)
    metric = 'diag'
    
    #Generate data from AR
    data = ar_.ar_data(A, er, nd)
    
    # Set parameters for analysis
    an_.set_params(nf = nf, ss = True, metric = metric,
                   detrend = True)
    
    #Call any connectivity routine routine. 
    #Here are some calling examples, uncomment your preferred one for use:
    
    # Set which connectivity method to use
    an_.set_params(alg = 'coh')
    # Call the method
    an_.measure_full(data)
    
    #Instead of setting the 'alg' parameter, you can call directly:
    #an_.coh_full(data), an_.pdc_full(data), etc.
        
    pp.show()
    
if __name__ == '__main__':
    
    simple_simulation_analysis()
    
    
