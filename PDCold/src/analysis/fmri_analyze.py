import time

from numpy import *
import scipy.stats as st
import scipy.signal as sig

from nifti import NiftiImage as ni

import analysis.explore as ex

import algorithms.pdc_alg as pdc

fil = sig.signaltools.lfilter
dr = sig.signaltools.detrend

print 'begin :', time.ctime()


nid = ni(r'D:\work\dados\fmri\resting\filter-reg\c1_+.ica\reg_standard\filtered_func_data.nii')
#nid = ni(r'D:\work\dados\fmri\resting\filter-reg\c2_setup7.ica\reg_standard\filtered_func_data.nii')
#nid = ni(r'D:\work\dados\fmri\resting\filter-reg\c2_setup7.ica\filtered_func_data.nii')
#nid = ni(r'D:\work\dados\fmri\resting\dados raw\c1')
#nid = ni(r'D:\work\dados\fmri\resting\ica\setup7.gica\groupmelodic.ica\melodic_IC')

mask = ni(r'D:\work\dados\fmri\resting\filter-reg\c1_+.ica\reg_standard\mask.nii').data.transpose()
#mask = mask[:,::-1,:]

data = nid.data

#data = data/std(data, axis=0)

data = data.transpose()

#datar = dr(data)

#lfb = array([1,1,1])/3.0
#lfa = [1.0]
#lfb = array([0.2])
#lfa = array([1.0, -0.8])
#lfb = array([0.5])
#lfa = array([1.0, -0.5])

#sh = data.shape

#dataf = fil(lfb, lfa, datar.reshape(-1,sh[-1])).reshape(sh)

#dataf[...,:2] = 0

#seed = dataf[71,14,6]
#seed = dataf[15,37,16]
#t6 = loadtxt(r'D:\work\dados\fmri\resting\ica\setup7.gica\groupmelodic.ica\report\t6.txt')
#seed = t6[:,1]

seed = data[19,9,20]
mask[19,9,20] = 0

func = lambda a, b: pdc.pdc([a, b], SS = False)
#func = st.pearsonr

#an = ex.seed_to_all(dataf, seed)
#an = ex.seed_to_all(data, seed, f = 'corr', mask = mask, default = -1)
an = ex.seed_to_all(data, seed, f = func, mask = mask)
an[19,9,20] = 1

#anni = ni(an[:,:,:,0].transpose())
#anni.save(r'D:\work\dados\fmri\resting\r2c1s19_9_20corr_mask.nii')

#anni = ni(an[:,:,:,1].transpose())
#anni.save(r'D:\work\dados\fmri\resting\r2c1s19_9_20pv_mask.nii')

anabs = abs(an)**2
anni = ni(anabs[:,:,:,0,1,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\r2c1s19_9_20pdc_0_1_all_mask.nii')

anni = ni(anabs[:,:,:,1,0,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\r2c1s19_9_20pdc_1_0_all_mask.nii')

print 'end: ', time.ctime()
