'''
Created on Jun 24, 2010

Simulate Integrate and fire network

@author: stein
'''

from pdc import *

from numpy import *

import scipy.stats as st
from scipy import randn, rand

def simIF(w, delay = 3, decay = 0.85, time = 100000, th = 1.0, noise_f = None):
    
    if (noise_f == None):
        stdn = 0.21
        #noise_f = lambda: st.norm.rvs(0, stdn)
        noise_f = lambda: rand()*stdn
    
    n = w.shape[0]
    
    v = zeros((n, time))
    s = zeros((n, time))
    for t in arange(delay, time):
        v[:,t] = v[:,t-1]*decay + dot(w,s[:,t-delay]) + noise_f()
        s[:,t] = v[:,t] > th
        v[v[:,t] > th, t] = 0
        
        if (t+1)%10000 == 0:
            print 'time', t+1
        
    return s, v

def test_net():
    n = 30
    
    w = 3*(1.0/n)*(randn(n,n) + 0.3)
    
    s, v = simIF(w)
    
    pr_.maxp = 30
    pr_.fixp = True
    
    s = pre_data(s)
    
    #A = ar_estim(s)
    A = 0
    
    print 'firerate', 1/(sum(s > 0, axis = 1)/100000.0)
    
    return A, s, v, w

def wtoa(w, delay = 3, decay = 0.85, p = 30):
    a = zeros(p)
    a[delay-1:] = w*decay**(arange(p-delay+1))
    return a

if __name__ == "__main__":
    
    test_net()