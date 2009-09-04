from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import scipy.stats as st
import time
from numpy.random import randn
from numpy.random import rand
from numpy.random import multivariate_normal as mnorm


import algorithms.asymp as ass_
import algorithms.pdc_alg as pdc_
from data_simulation.ar_data import ar_data
from data_simulation.ar_data import ar_models
from algorithms.ar_fit import nstrand

def gen_data_Ding(m):
    p = 3
    x = zeros([5,m])
    e = mnorm(zeros(7), identity(7), m).T
    b = 2*ones(5)
    c = 5*ones(5)
    a = rand(5)
    for i in arange(p, m):
        x[0,i] = 0.95*sqrt(2)*x[0,i-1] - 0.9025*x[0,i-2] + e[0,i] + a[0]*e[5,i] + b[0]*e[6,i-1] + c[0]*e[6,i-2]
        x[1,i] = 0.5*x[0,i-2] + e[1,i] + a[1]*e[5,i] + b[1]*e[6,i-1] + c[1]*e[6,i-2]
        x[2,i] = -0.4*x[0,i-3] + e[2,i]+ a[2]*e[5,i] + b[2]*e[6,i-1] + c[2]*e[6,i-2]
        x[3,i] = -0.5*x[0,i-2] + 0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[3,i]+ a[3]*e[5,i] + b[3]*e[6,i-1] + c[3]*e[6,i-2]
        x[4,i] = -0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[4,i]+ a[4]*e[5,i] + b[4]*e[6,i-1] + c[4]*e[6,i-2]

    return x

def teste_Ding():
    nd = 2000
    nf = 10
    alpha = 0.01
    n = 5
    maxp = 3
    metric = 'diag'
    
    #Generate data from AR
    data = gen_data_Ding(nd)
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = nstrand(data, maxp = maxp)
    #Calculate the connectivity and statistics
    mes, th, ic1, ic2 = ass_.asymp_pdc(data, Aest, nf, erest, 
                                   maxp, alpha = alpha, metric = metric)
    pdc_.plot_all(mes, th, ic1, ic2, nf = nf)

if __name__ == "__main__":
    teste_Ding()