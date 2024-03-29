# -*- coding:utf-8 -*-
"""
Created on 2009

@author: Carlos Stein

This file has examples of uses of the PDC library.

It has some examples reproducing related articles and some new ones.

It should be used to validate the library and to learn how to use it.

Please refer to the manual to further explanations on usage.
"""

from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import scipy.stats as st
import time
from numpy.random import randn
from numpy.random import rand
from numpy.random import multivariate_normal as mnorm
from scipy.integrate import odeint

from pdc import *

def teste_simples():
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
    nf = 40
    
    #error of the confidence interval and threshold
    alpha = 0.05
    
    #model order parameters
    n = A.shape[0]
    maxp = A.shape[2]
    
    #type of PDC used (refer to manual to see what it means)
    metric = 'gen'
    
    #Generate data from AR
    data = ar_data(A, er, nd)
    
    #Call any connectivity routine routine. 
    #Here are some calling examples, uncomment your preferred one for use:
    
    #pdc_.measure_and_plot(data, 'dtf', nf = nf, ss = True)
    #pdc_.pdc_and_plot(data, nf = nf, ss = True)
    #pdc_.pdc_full(data, nf = nf, ss = True, metric = metric)
    coh_full(data, nf = nf, ss = True, metric = metric,
                  detrend = True)
    
    pp.show()
    
    
    #If you want step by step, you can do it this way:
    
    #Estimate AR parameters with Nuttall-Strand
    #Aest, erest = ar_estim(data, maxp)
    
    #Calculate the connectivity and statistics
    #mes, th, ic1, ic2 = ass_.asymp_pdc(data, Aest, nf, erest, 
    #                               maxp, alpha = alpha, metric = metric)
    
    #Plot result
    #pdc_.plot_all(mes, th, ic1, ic2)
    
    
    #Another step-by-step way, without statistics:
    
    #Estimate AR parameters with Nuttall-Strand
    #Aest, erest = ar_estim(data, maxp)
    
    #Calculate the connectivity
    #mes = pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric)
    
    #Plot result
    #pdc_.pdc_plot(mes)
    
    

def gen_data_Guo(m, dummy = 100, a = None, bv = 2, cv = 5):
    '''Guo et. al, 2008. Uncovering interactions in the freq. domain'''
    p = 3
    x = zeros([5,m+dummy])
    e = mnorm(zeros(7), identity(7), m+dummy).T
    b = bv*ones(5)
    c = cv*ones(5)
    if a is None:
        a = rand(5)
        
#    print 'a', a
#    print 'b', b
#    print 'c', c
        
    for i in arange(p, m+dummy):
        x[0,i] = 0.95*sqrt(2)*x[0,i-1] - 0.9025*x[0,i-2] + e[0,i] + a[0]*e[5,i] + b[0]*e[6,i-1] + c[0]*e[6,i-2]
        x[1,i] = 0.5*x[0,i-2] + e[1,i] + a[1]*e[5,i] + b[1]*e[6,i-1] + c[1]*e[6,i-2]
        x[2,i] = -0.4*x[0,i-3] + e[2,i]+ a[2]*e[5,i] + b[2]*e[6,i-1] + c[2]*e[6,i-2]
        x[3,i] = -0.5*x[0,i-2] + 0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[3,i]+ a[3]*e[5,i] + b[3]*e[6,i-1] + c[3]*e[6,i-2]
        x[4,i] = -0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[4,i]+ a[4]*e[5,i] + b[4]*e[6,i-1] + c[4]*e[6,i-2]

    return x[:,dummy:]

def gen_data_Baccala(m, dummy = 100, a = None, bv = 2, cv = 5):
    '''Baccala and Sameshima, 2001. Biol. Cyber.'''
    gen_data_Guo(m, dummy, a = zeros(5), bv = 0, cv = 0)
    
def model_Eichler_2006():
    '''Eq. 26 Eichler 2006 on the DTF limiar'''
    
    a = zeros([2,2,4])
    a[0,0,:] = [7.0/8, -4.0/5, 4.0/5, -1.0/2]
    a[1,0,:2] = [10.0/20, 10.0/20]
    e = identity(2)
    
    return a,e 
    

def teste_Guo():
    nd = 2000
    nf = 8
    alpha = 0.01
    n = 5
    maxp = 3
    metric = 'gen'
    sample_f = 200
    n_boot = 100
    
    set_params(maxp = maxp, fixp = True, nf = nf, ss = False, sample_f = sample_f,
                   alpha = alpha, metric = metric, normalize = False)
    
#    data = gen_data_Guo(nd)
#    pdc_.pdc_full(data, do_plot = True)
#    pp.show()
#    return
    
    er = zeros([5, 5, nf])
    for i in arange(n_boot):
        #Generate data from Guo`s model
        data = gen_data_Guo(nd)
    
        mes, th, ic1, ic2 = pdc_full(data, do_plot = False)
        er += (mes > th)
        
        pr_.v = False
        
        if i%10 == 0:
            print i
        
    #mea = mean(mes, axis = 0)
    #smea = sort(mes, axis = 0)
    #th = smea[(1-alpha)*n_boot]
    
    #print mes[:,0,0,5], th[0,0,5]
    #er = sum(mes > th, axis = 0)/float(n_boot)
    er = er/float(n_boot)
    
    print mean(er, axis = -1)
    
    res_.mes = er
    #gl_.res_.th = th
    
    plot_all()
        
    #print pdc_.gct(data, maxp=3)
    
    pp.show()

def teste_sunspot_melanoma():
    #Generate data from AR @IndentOk
    y=array([[1936,  1.0, 0.9,  40],
        [ 1937, 0.8, 0.8, 115],
        [ 1938, 0.8, 0.8, 100],
        [ 1939, 1.4, 1.3,  80],
        [ 1940, 1.2, 1.4,  60],
        [ 1941, 1.0, 1.2,  40],
        [ 1942, 1.5, 1.7,  23],
        [ 1943, 1.9, 1.8,  10],
        [ 1944, 1.5, 1.6,  10],
        [ 1945, 1.5, 1.5,  25],
        [ 1946, 1.5, 1.5,  75],
        [ 1947, 1.6, 2.0, 145],
        [ 1948, 1.8, 2.5, 130],
        [ 1949, 2.8, 2.7, 130],
        [ 1950, 2.5, 2.9,  80],
        [ 1951, 2.5, 2.5,  65],
        [ 1952, 2.4, 3.1,  20],
        [ 1953, 2.1, 2.4,  10],
        [ 1954, 1.9, 2.2,   5],
        [ 1955, 2.4, 2.9,  10],
        [ 1956, 2.4, 2.5,  60],
        [ 1957, 2.6, 2.6, 190],
        [ 1958, 2.6, 3.2, 180],
        [ 1959, 4.4, 3.8, 175],
        [ 1960, 4.2, 4.2, 120],
        [ 1961, 3.8, 3.9,  50],
        [ 1962, 3.4, 3.7,  35],
        [ 1963, 3.6, 3.3,  20],
        [ 1964, 4.1, 3.7,  10],
        [ 1965, 3.7, 3.9,  15],
        [ 1966, 4.2, 4.1,  30],
        [ 1967, 4.1, 3.8,  60],
        [ 1968, 4.1, 4.7, 105],
        [ 1969, 4.0, 4.4, 105],
        [ 1970, 5.2, 4.8, 105],
        [ 1971, 5.3, 4.8,  80],
        [ 1972, 5.3, 4.8,  65]])
    data=y[:,[3,2]].transpose()
   
    nf = 64
    alpha = 0.01
    metric = 'orig'
    maxp = 10
   
    set_params(plot_labels = ['sun spots', 'melanome'])

   
    pdc_full(data, maxp = maxp, nf = nf, ss = True, 
                 alpha = alpha, metric = metric, normalize = True, 
                 stat = 'asymp')
    
    
    dtf_full(data, maxp = maxp, nf = nf, ss = True, 
                 alpha = alpha, metric = metric, normalize = True, 
                 stat = 'asymp')
    
    pp.show()


def gen_winterhalter_2005_van_der_Pol(n, dummy = 100, dt = 0.01):
    
    w = array([1.5,1.48,1.53,1.44])
    sg = 1.5
    mi = 5
    t = array([[0, 0.2, 0, 0],
               [0.2, 0, 0, 0.2],
               [0.2, 0, 0, 0.2],
               [0, 0.2, 0, 0]])
    
    data = zeros([4, n+dummy])
    x = zeros(4)
    x1 = zeros(4)
    x2 = zeros(4)
    for j in arange(1,n+dummy):
        n = randn(4)
        x = x+x1*dt
        x1 = x1+x2*dt
        x2 = mi*(1 - x**2)*x1 - w**2 * x + \
             sg*n/sqrt(dt) + dot(t, x) - sum(t,1)*x
        data[:,j] = x
        #x = xn
        #x1 = x1n
        #x2 = x2n
    
    data = data[:,dummy:]
    return data

#t = array([[0, 0.2, 0, 0],
#           [0.2, 0, 0, 0.2],
#           [0.2, 0, 0, 0.2],
#           [0, 0.2, 0, 0]])
#w = array([1.5,1.48,1.53,1.44])
#sg = 1.5
#
#mi = 5
#
#n = 0
#
#def odewinter_der(y, tm):
#
#    print 'ytm', y, tm
#    y, y1 = y.reshape(2,4)
#    
#    #nr = n[:,]
#    y2 = mi*(1 - y**2)*y1 - w**2 * y + \
#         sg*n + dot(t, y) - sum(t,1)*y
#    newy = concatenate((y1,y2)) 
#    
#    print 'new', newy
#    return newy 
#    
#def gen_winterhalter_2005_van_der_Pol_odeint(np, dummy = 100, dt = 0.01):
#    #esta com problemas por causa do random, eu acho.
#    n = randn(4,np+dummy)
#    data = odeint(odewinter_der, zeros(8), [0,1], mxstep = 10)#linspace(0,(np+dummy)*dt,np+dummy))
#    print data[:4, 100:140]
#    return data[:4,dummy:]

def teste_data():
    subs = 50
    nd = 50000*subs
    nf = 40
    alpha = 0.05
    #n = 5
    maxp = 10
    metric = 'orig'
    
    #Generate data from AR
    #data = gen_winterhalter_2005_van_der_Pol_odeint(nd, dt = 0.5/subs)
    data = loadtxt('D:/work/dados/simulation/pol50000_sub05_euler.txt')
    #data = subsample(data, subs)
    #data = subsample(data, 2)
       
    data = data[:,:4000]
    
    pdc_full(data, maxp = maxp, nf = nf, ss = True, 
                  metric = metric, alpha = alpha,
                  normalize = False, detrend = True, fixp = True, stat = 'asymp', n_boot = 100)
    
    #Estimate AR parameters with Nuttall-Strand
    #Aest, erest = ar_estim(data, maxp)
    #Calculate the connectivity and statistics
    #mes, th, ic1, ic2 = ass_.asymp_pdc(data, Aest, nf, erest, 
    #                               maxp, alpha = alpha, metric = metric)
    #pdc_.plot_all(mes, th, ic1, ic2, nf = nf)



if __name__ == "__main__":
    #artigo()
    #teste_simples()
    
    #teste_Guo()
    teste_sunspot_melanoma()
    #teste_data()
    #a = gen_winterhalter_2005_van_der_Pol(30, 30)
    #print a
    
