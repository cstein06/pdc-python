#from numpy import *
from scipy.io.matlab.mio import savemat
import os.path
import time
import pickle

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
        
        self.alg = 'pdc' 
        self.metric = 'diag'
        self.normalize = False  
        self.detrend = True   
        self.nf = 64   
        self.sample_f = 1   
        self.maxp = 30   
        self.fixp = False   
        self.ss = True   
        self.power = True   
        
        #Statistics specs
        
        self.alpha = 0.05   
        self.stat = 'asymp' 
        self.n_boot = 1000
        
        #Plotting specs
        
        self.do_plot = True
        self.plotf = None
        self.plot_diag = False
        self.plot_color = None
        self.plot_labels = None
        self.plot_title = None
        self.logss = True
        self.sqrtmes = False
        self.plot_th = True
        self.plot_ic = False
        self.plot_states = None
        self.state_colors = ['k', 'lightblue', 'darkblue', 'pink', 'red', 'green']
        
        #States specs
     
        self.window_size = None
        self.valid_states = None
        self.st_dict = None
        
        #Logging specs
        
        self.do_log = False
        self.log_matlab = False
        self.root_dir = 'C:/pdcpython/'
        self.log_string = 'current_log'
        self.mat_file = '%s%s'
        self.log_file = '%s%s.log'
        self.pic_file = '%s%s.pic'
        self.data_descr = 'Current log. Please give description.'  
        
        #States logging specs
        
        self.do_states_log = True
        self.stinput = 'current_state'
        self.stlog_pr = '%s%s_param.pic'
        self.stlog_res = '%s%s_res'
        self.stlog_mean = '%s%s_mean'
        
        #MISC
        
        self.time = None
        self.version = __pdc_version__
        self.v = True #verbose

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
    
    f = open(aux_pic, 'w')
    
    pickle.dump(pr_, f)    
    pickle.dump(res_, f)
    f.close()
    
    if pr_.log_matlab:
        savemat(aux_mat, {'result':res_.mes}, oned_as = 'row')
    
    
def load_results():
    global res_
    global pr_ 
    
    f = open(pr_.pic_file, 'r')
    
    pr_ = pickle.load(f)
    res_ = pickle.load(f)
    
    f.close()

    #return prn_, resn_

def log_params(file = None, **args):
    global pr_

    aux_log = pr_.log_file % (pr_.root_dir, pr_.log_string) 
    
    read_args(args)
    
    if file is not None:
        aux_log = file

    if pr_.v:
        print 'Logging the parameters to file:', aux_log  
    
    pr_.time = time.ctime()

    f = open(aux_log, 'w')
    
    pickle.dump(pr_, f)  
    
    f.close()
    
def load_params():
    global pr_ 
    
    f = open(pr_.pic_file, 'r')
    
    pr_ = pickle.load(f)
    
    f.close()

    #return prn_, resn_

def log_windows_results(stres, stmean, ststds, nstates, bind = False):
    
    pr_.stinput = pr_.stinput.replace('.txt', '') + '_' + pr_.alg
    
    if bind is True:
        pr_.stinput += '_bind'
    
    aux_res = pr_.stlog_res % (pr_.root_dir, pr_.stinput)
    aux_mean = pr_.stlog_mean % (pr_.root_dir, pr_.stinput)
    aux_pr = pr_.stlog_pr % (pr_.root_dir, pr_.stinput)
    
    print '\nLogging the raw results in file:', aux_res  
    print 'Logging the mean results in file:', aux_mean
    
    if os.path.isfile(aux_res + '.mat'):
        print '\nOverwriting results .mat file!'
            
    if os.path.isfile(aux_mean + '.mat'):
        print '\nOverwriting mean .mat file!'
    
    for i in range(len(stmean)):
        suf = ""
        if i > 0:
            suf = '_' + str(i)
        savemat(aux_res + suf,
                {'result':stres[i], 'time':time.ctime()}, oned_as = 'row')
    
        savemat(aux_mean + suf,
                {pr_.alg + '_mean':stmean[i], pr_.alg + '_stds':ststds[i], 
                'nstates':nstates,
                'time':time.ctime()}, oned_as = 'row')
        
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
