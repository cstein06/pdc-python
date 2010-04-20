# -*- coding:utf-8 -*-
"""
Created on 22/03/2010

@author: eu
"""

import pdc.sim_data as ard_
import pdc.analysis as pdc_
import pdc.examples as exa_
from pdc.params import *
import pdc.ar_fits as fit_
import pdc.asymp as asy_
import pdc.plotting as pl_

from numpy import *
from numpy.random import rand

import scipy.stats as st

import matplotlib.pyplot as pp

def sunspot():

    nd = 10000
    alpha = 0.01
    maxp = 5
    metric = 'diag'
    nboot = 5000
    nf = 64
    sample_f = 1
    
    A, er = ard_.ar_models(5)
    
    #data = ar_data(A, er, nd)
    data = ard_.ar_models(2)
    
    set_params(maxp = maxp, nf = nf, plot_labels = ['Sunspot', 'Melanoma'],
               logss = False, plot_ic = True, sample_f = sample_f,
               metric = metric, alpha = alpha, stat = 'asymp')
    
    res1 = pdc_.pdc_full(data)
    
    root = 'G:\\stein\\producao\\pdc congresso baccala\\'
    
    savetxt(root + 'ic1_1_2.txt', res1[2][0,1])
    savetxt(root + 'ic1_2_1.txt', res1[2][1,0])
    savetxt(root + 'ic2_1_2.txt', res1[3][0,1])
    savetxt(root + 'ic2_2_1.txt', res1[3][1,0])
    
    pp.show()
    return
    
    res2 = pdc_.pdc_full(data, stat = 'boot', n_boot = nboot)
    
    ratio = res1[1]/res2[1]
    
    res_.mes = ratio
    pl_.plot_all()
    
    pp.show()
    return
    
    res1 = pdc_.pdc_full(data)
    
    return
    pp.figure()
    
    
    res2 = pdc_.pdc_full(data, maxp = maxp, metric = metric, nf = nf,
                         alpha = alpha, stat = 'boot', n_boot = nboot)
    
    n=data.shape[0]
    
    pp.figure()
    x = arange(nf)/(2.0*nf)
    for i in arange(n):
        for j in arange(n):
            pp.subplot(n,n,i*n+j+1)
            if i == j: continue
            
            if i == 1 and j == 0: 
                pp.plot(x,res1[1][i,j])
                pp.plot(x,res2[1][i,j])
                continue
            pp.plot(x,res1[3][i,j]-res1[2][i,j])
            pp.plot(x,res2[3][i,j]-res2[2][i,j])
    
    pp.show()
    
    return res1, res2

def various_Guo(m = 5000):
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
    #a = zeros(5)
    print a
    
    #data = ard_.ar_data(A, er, nd)
    #data = loadtxt('D:\\work\\producao\\pdc congresso baccala\\ES57_09_02_09_medias_test.txt').T
    #data = data[:2, :5000]
    #data = ard_.ar_models(2)
    
    n = 5
    nd = 2000
    
    res = zeros([m, n, n, pr_.nf])
    
    for i in arange(m):
    
        data = exa_.gen_data_Guo(nd, a=a, bv = 0, cv = 0)
    
        res[i] = pdc_.measure(data)
        
        pr_.v = False
        if (i+1) % 20 == 0:
            print 'iter:', i+1
           
    pr_.do_plot = False
    
    pr_.v = True
    
    big = 1
    bins = 40
    
    
    data = exa_.gen_data_Guo(nd*big, a=a, bv = 0, cv = 0)
    #data2 = exa_.gen_data_Guo(nd, a=a)
    #A, er = fit_.nstrand(data, maxp = 3)
    #pdc, th, ic1, ic2 = asy_.asymp_pdc()
    pdc, th, ic1, ic2 = pdc_.measure_full(data)
    std = sqrt(big)*(pdc-ic1)/st.norm.ppf(1-pr_.alpha/2.0)    
        
    
    #pa = pdc[3,0,2]
    pa = mean(res[:,3,0,2])
    sa = std[3,0,2]
    
    fig = pp.figure()
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.norm.ppf(x, loc = pa, scale = sa)
    xmi = min(y.min(), res[:,3,0,2].min())
    xma = max(y.max(), res[:,3,0,2].max())
    
    ax = fig.add_subplot(121)
    
    pp.plot(sorted(res[:,3,0,2]), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'b')
    
    pp.xlabel('Estimated gPDC')
    pp.ylabel('Quantile for Normal')
    #pp.yticks()
    ax.yaxis.set_ticklabels(['0.001', '0.050', '0.250', '0.750', '0.950', '0.999'])
    ax.yaxis.set_ticks([y[0.001*m], y[0.050*m], y[0.250*m], y[0.750*m], y[0.950*m], y[0.999*m]])
    
    #fig = pp.figure()
    
    ax = fig.add_subplot(122)
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.chi2.ppf(x, pr_.patdf[4,0,2], scale = 1/(pr_.patden[4,0,2]*2*nd))
    xmi = min(y.min(), res[:,4,0,2].min())
    xma = max(y.max(), res[:,4,0,2].max())
    pp.plot(sorted(res[:,4,0,2]), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'b')
    
    
    pp.xlabel('Estimated gPDC')
    pp.ylabel('Quantile for weighted Chi-squares')
    #pp.yticks()
    tic = array(['0.001', '0.500', '0.750', '0.950', '0.999'])
    ax.yaxis.set_ticklabels(tic)
    ax.yaxis.set_ticks([y[0.001*m], y[0.500*m], y[0.750*m], y[0.950*m], y[0.999*m]])
    
    pp.figure()
    
    pp.hist(res[:,3,0,2], bins = bins, normed = True)
    
    
    #$x = linspace(pa-3*sa, pa+3*sa, 100)
    x = linspace(0.15, 0.35, 100)
    pp.plot(x, st.norm.pdf(x, loc = pa, scale = sa))
    
    pp.figure()
    
    pp.hist(res[:,4,0,2], bins = bins, normed = True)
    
    #std = sqrt(big)*(pdc-ic1)/st.norm.ppf(1-pr_.alpha/2.0)  
    std = big*th
    
    pa = pdc[4,0,2]
    sa = std[4,0,2]
    
    print 'naosig', pa, sa
    
    
#                patdf = sum(d)**2/sum(d**2)
#                patden = sum(d)/sum(d**2)
#                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
#                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                
    
    x = linspace(0, 0.01, 100)
    print 'pat', pr_.patdf[4,0,2]
    pp.plot(x, st.chi2.pdf(x, pr_.patdf[4,0,2], scale = 1/(pr_.patden[4,0,2]*2*nd)))
    
    #pp.plot(x, st.norm.pdf(x, loc = pa, scale = sa))
    
#    f = open('G:\\stein\\producao\\pdc congresso baccala\\py\\hist.pic', 'w')
#    
#    cPickle.dump(res[:,3:5,0,2], f)
#    cPickle.dump(pdc[3:5,0,2], f)
#    cPickle.dump(std[3,0,2], f)
#    cPickle.dump(pr_.patdf[4,0,2], f)
#    cPickle.dump(1/(pr_.patden[4,0,2]*2*nd), f)
#        
#    f.close()
#    
    pp.show()
    
    return res


def quant_Guo(m = 2000):
    '''Make histogram of a frequency for the gPDC in Guos model
    and compare with asymp
    '''
    
    #Definition of the MVAR model
    #A, er = ard_.ar_models(5)
    #nd = 2000
    
    pr_.alpha = 0.01
    
    pr_.maxp = 3
    pr_.fixp = True
    pr_.sample_f = 1
    pr_.ss = False
    
    pr_.plot_ic = True
    
    pr_.alg = 'pdc'
    pr_.metric = 'diag'
    
    pr_.nf = 5
    
    #a = rand(5)
    a = zeros(5)
    #a = array([ 0.59,  0.52,  0.72,  0.98,  0.66])
    print a
    
    #b = 2
    #c = 5
    b = 0
    c = 0
    
    #data = ard_.ar_data(A, er, nd)
    #data = loadtxt('D:\\work\\producao\\pdc congresso baccala\\ES57_09_02_09_medias_test.txt').T
    #data = data[:2, :5000]
    #data = ard_.ar_models(2)
    
    n = 5
    nd = 2000
    
    res = zeros([m, n, n, pr_.nf])
    
    for i in arange(m):
    
        data = exa_.gen_data_Guo(nd, a=a, bv = b, cv = c)
    
        res[i] = pdc_.measure(data)
        
        pr_.v = False
        if (i+1) % 20 == 0:
            print 'iter:', i+1
           
    pr_.do_plot = False
    
    pr_.v = True
    
    big = 100
    bins = 40
    
    print '\nnow big simulation for asymp pdc'
    
    data = exa_.gen_data_Guo(nd*big, a=a, bv = b, cv = c)
    #data2 = exa_.gen_data_Guo(nd, a=a)
    #A, er = fit_.nstrand(data, maxp = 3)
    #pdc, th, ic1, ic2 = asy_.asymp_pdc()
    pdc, th, ic1, ic2 = pdc_.measure_full(data)
    std = sqrt(big)*(pdc-ic1)/st.norm.ppf(1-pr_.alpha/2.0)    
        
    
    #pa = pdc[3,0,2]
    pa = mean(res[:,3,0,2])
    sa = std[3,0,2]
    
    fig = pp.figure()
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.norm.ppf(x, loc = pa, scale = sa)
#    xmi = min(y.min(), res[:,3,0,2].min())
#    xma = max(y.max(), res[:,3,0,2].max())
    xmi = res[:,3,0,2].min()
    xma = res[:,3,0,2].max()
    
    ax = fig.add_subplot(121)

    pp.plot(sorted(res[:,3,0,2]), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'b')
    
    pp.xlabel('Estimated gPDC')
    pp.ylabel('Quantile for Normal')
    #pp.yticks()
    ax.yaxis.set_ticklabels(['0.001', '0.050', '0.250', '0.750', '0.950', '0.999'])
    ax.yaxis.set_ticks([y[0.001*m], y[0.050*m], y[0.250*m], y[0.750*m], y[0.950*m], y[0.999*m]])
    
    #fig = pp.figure()
    
    ax = fig.add_subplot(122)
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.chi2.ppf(x, pr_.patdf[4,0,2], scale = 1/(pr_.patden[4,0,2]*2*nd))
    #xmi = min(y.min(), res[:,4,0,2].min())
    #xma = max(y.max(), res[:,4,0,2].max())
    xmi = res[:,4,0,2].min()
    xma = res[:,4,0,2].max()
    pp.plot(sorted(res[:,4,0,2]), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'b')
    
    
    pp.xlabel('Estimated gPDC')
    pp.ylabel('Quantile for weighted Chi-squares')
    #pp.yticks()
    tic = array(['0.001', '0.500', '0.750', '0.950', '0.999'])
    ax.yaxis.set_ticklabels(tic)
    ax.yaxis.set_ticks([y[0.001*m], y[0.500*m], y[0.750*m], y[0.950*m], y[0.999*m]])
#    
#    pp.figure()
    
    
    #pp.plot(x, st.norm.pdf(x, loc = pa, scale = sa))
    
#    f = open('G:\\stein\\producao\\pdc congresso baccala\\py\\hist.pic', 'w')
#    
#    cPickle.dump(res[:,3:5,0,2], f)
#    cPickle.dump(pdc[3:5,0,2], f)
#    cPickle.dump(std[3,0,2], f)
#    cPickle.dump(pr_.patdf[4,0,2], f)
#    cPickle.dump(1/(pr_.patden[4,0,2]*2*nd), f)
#        
#    f.close()
#    
    pp.show()
    
    return res

def Guo_error(m = 20):
    '''Calculate error on Guo model
    '''
    
    #Definition of the MVAR model
    #A, er = ard_.ar_models(5)
    #nd = 2000
    
    pr_.alpha = 0.01
    
    pr_.maxp = 3
    pr_.fixp = True
    pr_.sample_f = 1
    pr_.ss = False
    
    pr_.plot_ic = True
    
    pr_.alg = 'pdc'
    pr_.metric = 'diag'
    
    pr_.nf = 10
    
    #a = rand(5)
    #a = zeros(5)
    #print a
    
    #data = ard_.ar_data(A, er, nd)
    #data = loadtxt('D:\\work\\producao\\pdc congresso baccala\\ES57_09_02_09_medias_test.txt').T
    #data = data[:2, :5000]
    #data = ard_.ar_models(2)
    
    n = 5
    nd = 2000
    
    sumr = zeros([n, n, pr_.nf])
    sumo = zeros([n, n, pr_.nf])
    for i in arange(m):
    
        data = exa_.gen_data_Guo(nd)
    
        res = pdc_.measure_full(data, do_plot = False)
        
        sumr += res[0]>res[1]
        sumo += logical_or(res[0]>res[3], res[0]<res[2])
        
        pr_.v = False
        if (i+1) % 20 == 0:
            print 'iter:', i+1
           
    sumr /= float(m)
    sumo /= float(m)
           
    pr_.do_plot = False
    
    #print sumr, sumo
    
    pr_.v = True
    
    return sumr, sumo
    
    

def compare_Guo():
    '''compare result for different a,b,c.'''
    
    pr_.alpha = 0.01
    
    pr_.maxp = 3
    pr_.fixp = True
    pr_.sample_f = 1
    pr_.ss = False
    
    pr_.plot_ic = True
    
    pr_.alg = 'pdc'
    pr_.metric = 'diag'
    
    pr_.nf = 5
    
    #a = rand(5)
    a = zeros(5)
    #a = array([ 0.59,  0.52,  0.72,  0.98,  0.66])
    print a
    
    
    #b = 2
    #c = 5
    b = 0
    c = 0
    
    #data = ard_.ar_data(A, er, nd)
    #data = loadtxt('D:\\work\\producao\\pdc congresso baccala\\ES57_09_02_09_medias_test.txt').T
    #data = data[:2, :5000]
    #data = ard_.ar_models(2)
    
    n = 5
    nd = 2000
    
    m = 5
    
    res = zeros([m, n, n, pr_.nf])
    
    for i in arange(m):
    
        data = exa_.gen_data_Guo(nd, a=a, bv = b, cv = c)
        res[i] = pdc_.measure_full(data, do_plot = True)[0]
    
    pp.show()
    
    print res[:,3,0,2], mean(res[:,3,0,2])
    
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
    pass
    #sunspot()
    #histogram_Guo(2000)
    #Guo_error()
    quant_Guo()
    #compare_Guo()
    #winterhalter()
    


