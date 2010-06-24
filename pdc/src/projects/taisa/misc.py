'''
Created on May 24, 2010

@author: stein
'''
from numpy import *
from pdc import *

import projects.taisa.spikes as sp
import projects.data_io.tests as t
import scipy.signal as sig
import cPickle as cp

root = '/media/8c8a676c-a8cd-4a18-ae81-0ad35333149b/dados/dados taisa/'
spi = cp.load(open(root + 'spi.pic', 'rb'))
lf30 = cp.load(open(root + 'lf30.pic', 'rb'))
emg = cp.load(open(root + 'emg.pic', 'rb'))

p = sp.spiketotime(spi)

#r = t.test_plxio()

fd = 5
cut = 40
fs = 2000
nyq = fs/2
b,a = sig.filter_design.butter(fd,cut/float(nyq))
ssig = sig.lfilter(b,a,p)
fsig = sig.lfilter(b,a,lf30)

hb,ha = sig.filter_design.butter(fd,hcut/float(nyq), btype = 'high')

te = vstack(ssig, fsig)

subs = subsample(te)

pdc_and_plot(te, maxp = 30, fixp = True, ar_est = 'yw')

figure(4); c,d = csd(p, lf30, Fs = 2000, NFFT = 10000, detrend = detrend_mean)

