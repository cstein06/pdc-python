from numpy import *
from scipy.stats import pearsonr

from nifti import NiftiImage as ni

#pearson_ = lambda a,b: pearsonr(a,b)[1]

def seed_to_all_mask(nii, seed, f, mask, default):
 
    ns = nii.shape
    nii_ = nii.reshape(-1,ns[-1])
    n = nii_.shape[0]
    mask_ = mask.reshape(-1)
    
    
    #print f(seed, nii_[mask_.argmax()])
    sample = array(f(seed, nii_[mask_.argmax()]))
    rsize = sample.shape
    #print rsize
    resp = empty(concatenate(([n], rsize)), dtype=sample.dtype) 
    #for i in range(n):
    cont = 0
    i = 0
    while i < n: # and cont < 100:
        if mask_[i]:
            resp[i] = f(seed, nii_[i])
            cont = cont+1
            if (cont%1000 == 0):
                print cont
        else:
            resp[i] = default
        i = i+1
    
    return resp.reshape(concatenate((ns[:-1],rsize)))
    

def seed_to_all(nii, seed, f = "corr", mask = None, default = 0.0):
    
    if f == "corr":
        f = lambda a,b: pearsonr(a,b)[0]
    if f == "pv":
        f = lambda a,b: pearsonr(a,b)[1]
    if f == "cp":
        f = lambda a,b: array(pearsonr(a,b))
        
    
    if mask != None:
        return seed_to_all_mask(nii, seed, f, mask, default)
    
    ns = nii.shape
    nii_ = nii.reshape(-1,ns[-1])
    n = nii_.shape[0]
    
    sample = f(seed, nii_[0])
    rsize = sample.shape
    resp = empty(concatenate(([n],rsize)), dtype = sample.dtype) 
    for i in range(n):
        resp[i] = f(seed, nii_[i])
    
    return resp.reshape(concatenate((ns[:-1],rsize)))
    
    
def seed_cube(nii, ia, ib, ic):
    s = zeros(nii.shape[3])
    for i in range(ia-1, ia+2):
        for j in range(ib-1, ib+2):
            for k in range(ic-1, ic+2):
                s = s+nii[i,j,k]
    return s/27.

def z_subsample():
#    for i in arange(1,7):
#        da = ni(r'D:\work\dados\fmri\resting\new\c' + str(i) + r'_new.ica\reg_standard\filtered_func_data.nii').data
#        nz = da.shape[1]
#        ind = int32(linspace(0,nz-1,nz/3))
#        da = da[:,ind]
#        da = ni(da)
#        da.save(r'D:\work\dados\fmri\resting\new\c' + str(i) + r'_new.ica\reg_standard\filtered_func_data_z3.nii')
#        da = ni(r'D:\work\dados\fmri\resting\new\c' + str(i) + r'_new.ica\reg_standard\mask.nii').data
#        nz = da.shape[0]
#        ind = int32(linspace(0,nz-1,nz/3))
#        da = da[ind]
#        da = ni(da)
#        da.save(r'D:\work\dados\fmri\resting\new\c' + str(i) + r'_new.ica\reg_standard\mask_z3.nii')
#    
#    da = ni(r'D:\work\dados\fmri\resting\new\new.gica\groupmelodic.ica\stats\thresh_zstat4.nii').data
#    nz = da.shape[0]
#    ind = int32(linspace(0,nz-1,nz/3))
#    da = da[ind]
#    da = ni(da)
#    da.save(r'D:\work\dados\fmri\resting\new\new.gica\groupmelodic.ica\stats\thresh_zstat4_z3.nii')
    
    da = ni(r'D:\work\dados\fmri\resting\new\new.gica\groupmelodic.ica\stats\thresh_zstat6.nii').data
    nz = da.shape[0]
    ind = int32(linspace(0,nz-1,nz/3))
    da = da[ind]
    da = ni(da)
    da.save(r'D:\work\dados\fmri\resting\new\new.gica\groupmelodic.ica\stats\thresh_zstat6_z3.nii')
  
def calc_mean():
    zz,aa,bb,cc = 64,15,54,45
    sumpdc01 = zeros([zz,aa,bb,cc])
    sumpdc10 = zeros([zz,aa,bb,cc])
    sumss0 = zeros([zz,aa,bb,cc])
    sumss1 = zeros([zz,aa,bb,cc])
    sumcoh = zeros([zz,aa,bb,cc])
    
    for i in arange(1,7):
       
        a = ni(r'D:\work\dados\fmri\resting\results\newc' + str(i) + r's13_25_11pdc_0_1_z3.nii')
        sumpdc01 = sumpdc01 + a.data
        del a        
        
        a = ni(r'D:\work\dados\fmri\resting\results\newc' + str(i) + r's13_25_11pdc_1_0_z3.nii')
        sumpdc10 = sumpdc10 + a.data
        del a        
        
        a = ni(r'D:\work\dados\fmri\resting\results\newc' + str(i) + r's13_25_11ss_0_z3.nii')
        sumss0 = sumss0 + a.data
        del a        
        
        a = ni(r'D:\work\dados\fmri\resting\results\newc' + str(i) + r's13_25_11ss_1_z3.nii')
        sumss1 = sumss1 + a.data
        del a        
        
        a = ni(r'D:\work\dados\fmri\resting\results\newc' + str(i) + r's13_25_11coh_z3.nii')
        sumcoh = sumcoh + a.data
        del a        
        
    a = ni(sumpdc01/6)
    a.save(r'D:\work\dados\fmri\resting\results\newmeans13_25_11pdc_0_1_z3.nii')
    del a
    
    a = ni(sumpdc10/6)
    a.save(r'D:\work\dados\fmri\resting\results\newmeans13_25_11pdc_1_0_z3.nii')
    del a
    
    a = ni(sumss0/6)
    a.save(r'D:\work\dados\fmri\resting\results\newmeans13_25_11ss_0_z3.nii')
    del a
    
    a = ni(sumss1/6)
    a.save(r'D:\work\dados\fmri\resting\results\newmeans13_25_11ss_1_z3.nii')
    del a
    
    a = ni(sumcoh/6)
    a.save(r'D:\work\dados\fmri\resting\results\newmeans13_25_11coh_z3.nii')
    del a

if __name__ == "__main__":
    #z_subsample()
    calc_mean()
    pass
    
      