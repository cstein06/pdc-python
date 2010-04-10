from numpy import *
from scipy.stats import pearsonr

#pearson_ = lambda a,b: pearsonr(a,b)[1]

def seed_to_all_mask(nii, seed, f, mask, default):
 
    ns = nii.shape
    nii_ = nii.reshape(-1,ns[-1])
    n = nii_.shape[0]
    mask_ = mask.reshape(-1)
    
    #print f(seed, nii_[mask_.argmax()])
    sample = f(seed, nii_[mask_.argmax()])
    rsize = sample.shape
    #print rsize
    resp = empty(concatenate(([n], rsize)), dtype=sample.dtype) 
    #for i in range(n):
    cont = 0
    i = 0
    while i < n:
        if mask_[i]:
            resp[i] = f(seed, nii_[i])
            cont = cont+1
            if (cont%1000 == 0):
                print cont
        else:
            resp[i] = ones(rsize)*default
        i = i+1
    
    return resp.reshape(concatenate((ns[:-1],rsize)))
    

def seed_to_all(nii, seed, f = "corr", mask = None, default = 0.0):
    
    if f == "corr":
        f = lambda a,b: pearsonr(a,b)[0]
    if f == "pv":
        f = lambda a,b: pearsonr(a,b)[1]
    
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
    
    