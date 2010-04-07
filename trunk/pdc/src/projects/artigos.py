# -*- coding:utf-8 -*-
"""
Created on 22/03/2010

@author: eu
"""

import pdc.ar_data as ard_
import pdc.analysis as pdc_
import pdc.examples as exa_
from pdc.globals import *
import pdc.ar_fit as fit_
import pdc.asymp as asy_

from numpy import *
from numpy.random import rand

import scipy.stats as st

import matplotlib.pyplot as pp

def histogram_Guo(m = 5000):
    '''Make histogram of a frequency for the gPDC in Guos model
    and compare with asymp
    '''
    
    #Definition of the MVAR model
    #A, er = ard_.ar_models(5)
    #nd = 2000
    
    pr_.alpha = 0.05
    
    pr_.maxp = 3
    pr_.fixp = True
    pr_.sample_f = 1
    pr_.ss = False
    
    pr_.plot_ic = True
    
    pr_.alg = 'pdc'
    pr_.metric = 'diag'
    
    pr_.nf = 5
    
    a = rand(5)
    print a
    
    #data = ard_.ar_data(A, er, nd)
    #data = loadtxt('D:\\work\\producao\\pdc congresso baccala\\ES57_09_02_09_medias_test.txt').T
    #data = data[:2, :5000]
    #data = ard_.ar_models(2)
    
    n = 5
    nd = 2000
    
    res = zeros([m, n, n, pr_.nf])
    
    for i in arange(m):
    
        data = exa_.gen_data_Guo(nd, a=a)
    
        res[i] = pdc_.measure(data)
        
        pr_.v = False
        if (i+1) % 20 == 0:
            print 'iter:', i+1
           
    pr_.do_plot = False
    
    pr_.v = True
    
    big = 400
    bins = 40
    
    
    data = exa_.gen_data_Guo(nd*big, a=a)
    #data2 = exa_.gen_data_Guo(nd, a=a)
    #A, er = fit_.nstrand(data, maxp = 3)
    #pdc, th, ic1, ic2 = asy_.asymp_pdc()
    pdc, th, ic1, ic2 = pdc_.measure_full(data)
    std = sqrt(big)*(pdc-ic1)/st.norm.ppf(1-pr_.alpha/2.0)    
        
    pp.hist(res[:,3,0,2], bins = bins, normed = True)
    
    pa = pdc[3,0,2]
    sa = std[3,0,2]
    #$x = linspace(pa-3*sa, pa+3*sa, 100)
    x = linspace(0, 0.6, 100)
    pp.plot(x, st.norm.pdf(x, loc = pa, scale = sa))
    
    pp.figure()
    
    pp.hist(res[:,4,0,2], bins = bins, normed = True)
    
    #std = sqrt(big)*(pdc-ic1)/st.norm.ppf(1-pr_.alpha/2.0)  
    std = big*th
    
    pa = pdc[4,0,2]
    #sa = std[4,0,2]
    
    
#                patdf = sum(d)**2/sum(d**2)
#                patden = sum(d)/sum(d**2)
#                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
#                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                
    
    x = linspace(0, 0.1, 100)
    pp.plot(x, st.chi2.pdf(x, pr_.patdf[4,0,2], scale = 1/(pr_.patden[4,0,2]*2*nd)))
    
    f = open('G:\\stein\\producao\\pdc congresso baccala\\py\\hist.pic', 'w')
    
    pickle.dump(res[:,3:5,0,2], f)
    pickle.dump(pdc[3:5,0,2], f)
    pickle.dump(std[3,0,2], f)
    pickle.dump(pr_.patdf[4,0,2], f)
    pickle.dump(1/(pr_.patden[4,0,2]*2*nd), f)
        
    f.close()
    
    pp.show()
    
    return res
    
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
    
    histogram_Guo()
    #winterhalter()
    
