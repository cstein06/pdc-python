import time

from numpy import *
import scipy.stats as st
import scipy.signal as sig

from nifti import NiftiImage as ni
import analysis.explore as ex
import algorithms.pdc_alg as pdc

print 'begin :', time.ctime()

data = ni(r'D:\work\dados\fmri\resting\new\c2_new.ica\reg_standard\filtered_func_data_z3.nii').data.transpose()
mask = ni(r'D:\work\dados\fmri\resting\new\c2_new.ica\reg_standard\mask_z3.nii').data.transpose()

aa,bb,cc = 13,25,11
seed = ex.seed_cube(data,aa,bb,cc)
mask[aa,bb,cc] = 0

#func = lambda a, b: pdc.coh([a, b])
func = lambda a, b: pdc.pdc_ss_coh([a, b])
#func = lambda a, b: pdc.pdc([a, b], SS = False)



#an = ex.seed_to_all(dataf, seed)
#an = ex.seed_to_all(data, seed, f = 'corr', mask = mask, default = -1)
#an = ex.seed_to_all(data, seed, f = 'cp', mask = mask, default = array([0,1]))
an = ex.seed_to_all(data, seed, f = func, mask = mask, default = 0)
an[aa,bb,cc] = 2

print an.shape
#an = abs(an[:23])
#an = an**2

#anni = ni(an[:,:,:,0].transpose())
#anni.save(r'D:\work\dados\fmri\resting\results\newc1s19_9_20corr.nii')

#anni = ni(an[:,:,:,1].transpose())
#anni.save(r'D:\work\dados\fmri\resting\results\newc1s19_9_20pv.nii')

anni = ni(an[:,:,:,0,0,1,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\results\newc2s13_25_11pdc_0_1_z3.nii')

anni = ni(an[:,:,:,0,1,0,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\results\newc2s13_25_11pdc_1_0_z3.nii')

anni = ni(an[:,:,:,1,0,0,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\results\newc2s13_25_11ss_0_z3.nii')

anni = ni(an[:,:,:,1,1,1,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\results\newc2s13_25_11ss_1_z3.nii')

anni = ni(an[:,:,:,2,0,1,:].transpose())
anni.save(r'D:\work\dados\fmri\resting\results\newc2s13_25_11coh_z3.nii')


print 'end: ', time.ctime()
