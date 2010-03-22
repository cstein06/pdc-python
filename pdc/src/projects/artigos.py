# -*- coding:utf-8 -*-
"""
Created on 22/03/2010

@author: eu
"""

import pdc.ar_data as ard_
import pdc.analysis as pdc_
import pdc.examples as exa_

def teste_simples():
    '''Simple test of connectivity routines
    
    Here we simulate some data using the ar_data module, from a given
    MVAR matrix A and covariance matrix er.
    
    '''
    
    #Definition of the MVAR model
    A, er = ard_.ar_models(5)
    
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
    metric = 'euc'
    
    #Generate data from AR
    data = ard_.ar_data(A, er, nd)
    
    #Call any connectivity routine routine. 
    #Here are some calling examples, uncomment your preferred one for use:
    
    #pdc_.measure_and_plot(data, 'dtf', nf = nf, ss = True)
    #pdc_.pdc_and_plot(data, nf = nf, ss = True)
    pdc_.pdc_full(data, nf = nf, ss = False, metric = metric)
    #pdc_.coh_full(data, nf = nf, ss = True, metric = metric,
    #              detrend = True)
    
    
    #If you want step by step, you can do it this way:
    
    #Estimate AR parameters with Nuttall-Strand
    #Aest, erest = ar_fit(data, maxp)
    
    #Calculate the connectivity and statistics
    #mes, th, ic1, ic2 = ass_.asymp_pdc(data, Aest, nf, erest, 
    #                               maxp, alpha = alpha, metric = metric)
    
    #Plot result
    #pdc_.plot_all(mes, th, ic1, ic2)
    
    
    #Another step-by-step way, without statistics:
    
    #Estimate AR parameters with Nuttall-Strand
    #Aest, erest = ar_fit(data, maxp)
    
    #Calculate the connectivity
    #mes = pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric)
    
    #Plot result
    #pdc_.pdc_plot(mes)
    
def winterhalter():
    
    subs = 50
    
    nd = 10000*subs
    
    metric = 'diag'
    
    nf = 64
    
    data = exa_.gen_winterhalter_2005_van_der_Pol(nd)
    
    data = data[:,::subs]
    
    print data.shape
    
    pdc_.pdc_full(data, nf = nf, ss = False, metric = metric)
    
if __name__ == "__main__":
    
    #teste_simples()
    winterhalter()
    
