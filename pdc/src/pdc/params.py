__all__ = ['pr_', 'res_', 'reset', 'log_results', 'load_results',
           'load_params', 'load_params', 'log_windows_results',
           'read_args', 'read_results', 'set_params', 'set_results']

#from numpy import *
from scipy.io.matlab.mio import savemat
import os.path
import time
import cPickle

__pdc_version__ = 0.1

mnames_ = {'coh':'Coherence',
          'pdc':'Partial Directed Coherence',
          'dtf':'Directed Transfer Function',
          'ss':'Spectral Density',
          'pc': 'Partial Coherence',
          'ar': 'VAR model'}


class Param():  
    def __init__(self):
        
        #Analysis specs
        
        self.alg = 'pdc' # Analysis method
        self.metric = 'diag' # PDC normalization
        self.normalize = False # Normalize data to var = 1
        self.detrend = True # Remove linear trend from data
        self.nf = 64 # Number of frequency points
        self.sample_f = 1 # Data sampling frequency
        self.maxp = 30 # Max order of estimated VAR
        self.fixp = False # Fix VAR order at maxp or not
        self.test_allp = False # Find AIC global minimum VAR order, not first minimum
        self.ss = True # Calculate power spectrum also
        self.power = True # Calculate power of the measure (instead of complex value)
        self.ar_estim = 'ns' # VAR estimator: ns (Nutall-Strand) or yw (Yule-Walker)
        self.reuse_A = False # Use the res_.A and res_.er, don't do new ar_estim
        
        #Statistics specs
        
        self.alpha = 0.05 # Alpha error for statistics 
        self.stat = 'asymp' # Bootstrap (boot) or Asymptotic (asymp) statistics
        self.n_boot = 1000 # Number of bootstrap repetitions
        
        #Plotting specs
        
        self.do_plot = True # Plot result
        self.plotf = None # Plot until this frequency
        self.plot_diag = False # Plot diagonal results (lots of times not needed)
        self.plot_color = None # Plot color to use
        self.plot_labels = None # Labels for each signal
        self.plot_title = None # Plot title
        self.logss = True # Plot log of power spectrum
        self.sqrtmes = False # Take sqrt of results
        self.plot_th = True # Plot threshold
        self.plot_ic = False # Plot confidence intervals
        
        #Logging specs
        
        self.do_log = False # Log results
        self.log_matlab = False # Log matlab format
        self.root_dir = 'G:\\stein\\dados\\edu_comp\\' # root data directory
        self.log_string = 'current_log' # file string
        self.mat_file = '%s%s' # matlab file name
        self.log_file = '%s%s.log' # log file name
        self.pic_file = '%s%s.pic' # cPickle file name
        self.data_descr = 'Current log. Please give description.' # data description
        
        # States plotting specs
        
        self.plot_states = None # Which states to plot
        self.state_colors = ['k', 'lightblue', 'darkblue', 
                             'pink', 'red', 'green'] # Plot color of each state
        
        #States specs
     
        self.window_size = 1 # Size of data windows
        self.valid_states = [1,2,3,4,5,6] # States to be analysed
        self.st_dict = None # Do not edit this
        
        #States logging specs
        
        self.do_states_log = True # Log results 
        self.output_dir = 'G:\\stein\\dados\\edu_comp\\results\\' # Output dir
        self.stinput = 'current_state' # file string
        self.stlog_pr = '%s%s_param.pic' # params cPickle string
        self.stlog_res = '%s%s_res' # matlab result string
        self.stlog_mean = '%s%s_mean' # matlab mean result string
        self.stpic_file = '%s%s_res.pic' # cPickle result string
        
        #MISC
        
        self.v = True # verbose. False for no text output.
        self.time = None
        self.version = __pdc_version__

pr_ = Param()

class Results():  
    def __init__(self):
        self.data_shape = None
        self.A = None  
        self.er = None   
        self.mes = None   
        self.ss = None 
        self.alg = None   
        self.metric = None 
        self.p = None   
        self.Af = None 
        self.th = None
        self.ic1 = None   
        self.ic2 = None
            

res_ = Results()


def reset():
    global res_
    global pr_
    
    res_ = Results()
    pr_ = Param()

#
#pm_ = {}
#pm_['metric'] = 'diag'
#pm_['normalize'] = False  
#pm_['detrend'] = True   
#pm_['nf'] = 64   
#pm_['sample_f'] = 1   
#pm_['maxp'] = 30   
#pm_['fixp'] = False   
#pm_['alg'] = 'pdc' 
#pm_['logss'] = True
#pm_['ss'] = True   
#pm_['power'] = True   
#pm_['alpha'] = 0.05   
#pm_['stat'] = 'asymp' 
#pm_['n_boot'] = 1000
#pm_['plotf'] = None
#
#metric_ = 'diag'
#normalize_ = False  
#detrend_ = True   
#nf_ = 64   
#sample_f_ = 1   
#maxp_ = 30   
#fixp_ = False   
#alg_ = 'pdc' 
#logss_ = True
#ss_ = True   
#power_ = True   
#alpha_ = 0.05   
#stat_ = 'asymp' 
#n_boot_ = 1000
#plotf_ = None


def log_results(**args):
    
    global res_
    global pr_

    aux_log = pr_.log_file % (pr_.root_dir, pr_.log_string) 
    aux_pic = pr_.pic_file % (pr_.root_dir, pr_.log_string) 
    aux_mat = pr_.mat_file % (pr_.root_dir, pr_.log_string) 
    
    read_args(args)

    if pr_.v:
        print 'Logging the results in file:', aux_log  
    
    pr_.time = time.ctime()

    f = open(aux_log, 'w')
    
    f.write(pr_.data_descr + '\n\n')
    
    f.write('Parameters:' + '\n\n')
    f.write(str(vars(pr_)) + '\n\n')
    
    f.write('Results:' + '\n\n')
    f.write(str(vars(res_)) + '\n\n')
    
    f.close()
    
    f = open(aux_pic, 'wb')
    
    cPickle.dump(pr_, f)    
    cPickle.dump(res_, f)
    f.close()
    
    if pr_.log_matlab:
        savemat(aux_mat, {'result':res_.mes}, oned_as = 'row')
    
    
def load_results():
    global res_
    global pr_ 
    
    f = open(pr_.pic_file, 'rb')
    
    pr_ = cPickle.load(f)
    res_ = cPickle.load(f)
    
    f.close()

    #return prn_, resn_

def log_params(file = None, **args):
    global pr_

    aux_log = pr_.log_file % (pr_.root_dir, pr_.log_string) 
    
    read_args(args)
    
    if file is not None:
        aux_log = file

    if pr_.v:
        print '\nLogging the parameters to file:', aux_log  
    
    pr_.time = time.ctime()

    f = open(aux_log, 'wb')
    
    cPickle.dump(pr_, f)  
    
    f.close()
    
def load_params(file):
    global pr_ 
    
    f = open(pr_.output_dir + file, 'rb')
    
    pr_ = cPickle.load(f)
    
    f.close()

    #return prn_, resn_

#def load_win_pic(file):
#    
#    print pr_.output_dir + file
#    f = open(pr_.output_dir + file, 'r')
#    
#    stres = cPickle.load(f)    
#    #stmean = cPickle.load(f)  
#    #ststds = cPickle.load(f)  
#    #nstates = cPickle.load(f)
#    ststds = None  
#    stmean = None
#    nstates= None
#    
#    #print stres.shape
#    #print stmean.shape
#    #print ststds.shape
#        
#    f.close()
#    
#    return stres, stmean, ststds, nstates

def log_windows_results(stres, stmean, ststds, nstates, states, bind = False):
    
    if pr_.output_dir is None:
        aux_dir = pr_.root_dir
    else:
        aux_dir = pr_.output_dir
    
    stiaux = pr_.stinput.replace('.txt', '') + '_' + pr_.alg
    
    if bind is True:
        pr_.stinput += '_bind'
    
    aux_res = pr_.stlog_res % (aux_dir, stiaux)
    aux_mean = pr_.stlog_mean % (aux_dir, stiaux)
    aux_pr = pr_.stlog_pr % (aux_dir, stiaux) 
    aux_pic = pr_.stpic_file % (aux_dir, stiaux)
    
    print '\nLogging the raw results in file:', aux_res  
    print 'Logging the mean results in file:', aux_mean
    print 'Logging the cPickled results in file:', aux_pic
    
    if os.path.isfile(aux_res + '.mat'):
        print '\nOverwriting results .mat file!'
            
    if os.path.isfile(aux_mean + '.mat'):
        print '\nOverwriting mean .mat file!'
            
    if os.path.isfile(aux_pic):
        print '\nOverwriting pic file!'
    
    for i in range(len(stmean)):
        suf = ""
        if i > 0:
            suf = '_' + str(i)
        savemat(aux_res + suf,
                {'result':stres[i], 'time':time.ctime()}, oned_as = 'row')
    
        savemat(aux_mean + suf,
                {pr_.alg + '_mean':stmean[i], pr_.alg + '_stds':ststds[i], 
                'nstates':nstates, 'states':states,
                'time':time.ctime()}, oned_as = 'row')
        
    f = open(aux_pic, 'wb')
    
    cPickle.dump(stres[0], f)    
    cPickle.dump(stmean[0], f)  
    cPickle.dump(ststds[0], f)  
    cPickle.dump(nstates, f)
    cPickle.dump(states, f)
        
    f.close()
    
        
    log_params(file = aux_pr)
    
def read_args(args):
    global res_
    global pr_ 
    
    for a,b in args.items():
        if not vars(pr_).has_key(a):
            print "Parameter given misspelled:", a
        vars(pr_)[a] = b

        #exec 'pr_.' + a + ' = b'
        
def read_results(args):
    
    for a,b in args.items():
        if not vars(res_).has_key(a):
            print "Parameter given misspelled"
        vars(res_)[a] = b
        #exec 'res_.' + a + ' = b'

def set_params(**args):
    read_args(args)
    
def set_results(**args):
    read_results(args)
