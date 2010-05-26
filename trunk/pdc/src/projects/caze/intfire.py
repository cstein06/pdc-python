'''
Created on Apr 28, 2010

@author: Marcelo Gomes Mattar
'''


from numpy import *
import scipy.stats as st
from numpy.matlib import randn
from numpy.random import rand
import matplotlib.pyplot as pp
import cProfile
import pdc.analysis as an


def bleh(n, alpha, beta, simlength):
    t = 3
    v = zeros((n,simlength),float)
    u = zeros((n,simlength),float)
    A = randn(n,n)/10
#    A = (1/n)*ones((n,n),float)
#    v[:,t-1] = 0.9*ones((n,1),float)
    while t < simlength-1:          
        v[:,t] = alpha*v[:,t-1] + dot(A,u[:,t-3]) + st.bernoulli.rvs(beta,size=n)
        for i in range(0,n):   
            if v[i,t-1] >= 1:
                u[i,t] = 1
                v[i,t] = -0.3
            else:
                u[i,t] = 0
            
        t = t + 1
            
#    return u[0,:]
    return u[0,:]



def network(W, Vrest = -65, Vthr = -50, Vreset = -70, noisestd = 30, spiketocurrent = 300, Rin = 10, tau = 15, timestep = 0.1, simtime = 100):
    nsteps = simtime/timestep
    n = W.shape[1]
    Vmem = Vrest * ones([n,nsteps])
    prevVmem = Vrest * ones(n)
    spiketrain = zeros([n,nsteps])
    Iin  = zeros(n)

    for t in arange(1,nsteps):
        Iin = spiketocurrent*dot(spiketrain[:,t-1],W) + noisestd*randn(n)
        Vmem[:,t] = prevVmem + ((Vrest - prevVmem) + Rin*Iin)*(timestep/tau)
        for i in arange(0,n):
            if Vmem[i,t] > Vthr:
                if t < (5/timestep):
                    if sum(spiketrain[i,0:(t-1)]) >= 1:
                        Vmem[:,t] = prevVmem + (Vrest - prevVmem)*(timestep/tau)   
                        prevVmem[i] = Vmem[i,t]
                        spiketrain[i,t] = 0
                    else:
                        prevVmem[i] = Vreset
                        spiketrain[i,t] = 1
                        Vmem[i,t] = Vthr
                else:
                    if sum(spiketrain[i,(t-(5/timestep)):(t-1)]) >= 1:
                        Vmem[:,t] = prevVmem + (Vrest - prevVmem)*(timestep/tau)   
                        prevVmem[i] = Vmem[i,t]
                        spiketrain[i,t] = 0
                    else:
                        prevVmem[i] = Vreset
                        spiketrain[i,t] = 1
                        Vmem[i,t] = Vthr
            else:
                prevVmem[i] = Vmem[i,t]
                spiketrain[i,t] = 0
            
    return Vmem, spiketrain


if __name__ == '__main__':
#    simul = neuron(100, 0.9, 0.01, 1000)
#    pp.plot(simul.T)
#    pp.show()
    n = 3
    W = zeros([n,n])
    W[0,1] = 1
    W[2,1] = 1
#    W = rand(3,3)
    W[arange(n),arange(n)] = 0
    print W
    Vmem,spiketrain = network(W)
    
#    print spiketrain
#    print Vmem
#    res = an.pdc_full([Vmem[0,:],Vmem[1,:],Vmem[2,:]])
    res = an.pdc_full([spiketrain[0,:],spiketrain[1,:],spiketrain[2,:]])
    pp.show()
    
    #===========================================================================
    # pp.figure(1)
    # pp.subplot(231)
    # pp.plot(Vmem[0,:].T)
    # pp.subplot(232)
    # pp.plot(Vmem[1,:].T)
    # pp.subplot(233)
    # pp.plot(Vmem[2,:].T)
    # pp.subplot(234)
    # pp.plot(spiketrain[0,:].T)
    # pp.subplot(235)
    # pp.plot(spiketrain[1,:].T)
    # pp.subplot(236)
    # pp.plot(spiketrain[2,:].T)
    # pp.show()
    #===========================================================================
    