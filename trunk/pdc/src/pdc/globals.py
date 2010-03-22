#from numpy import *
import pickle

mnames_ = {'coh':'Coherence',
          'pdc':'Partial Directed Coherence',
          'dtf':'Directed Transfer Function',
          'ss':'Spectral Density',
          'pc': 'Partial Coherence'}

state_colors = ['k', 'lightblue', 'darkblue', 'pink', 'red', 'green']

class Param():  
    def __init__(self):
        self.metric = 'diag'
        self.normalize = False  
        self.detrend = True   
        self.nf = 64   
        self.sample_f = 1   
        self.maxp = 30   
        self.fixp = False   
        self.alg = 'pdc' 
        self.logss = False
        self.sqrtmes = False
        self.ss = True   
        self.power = True   
        self.alpha = 0.05   
        self.stat = 'asymp' 
        self.n_boot = 1000
        self.plotf = None
        self.window_size = None
        self.plot_color = None
        self.do_plot = True
        self.do_log = True
        self.do_window_log = True
        self.v = True #verbose

pr_ = Param()

class Results():  
    def __init__(self):
        self.log_file = 'current_log.log'
        self.pic_file = 'current_log.pic'
        self.data_descr = 'Current log. Please give description.'  
        self.data_shape = None
        self.A = None  
        self.er = None   
        self.mes = None   
        self.ss = None 
        self.alg = None   
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

    read_args(args)

    f = open(res_.log_file, 'w')
    
    f.write(res_.data_descr + '\n\n')
    
    f.write('Parameters:' + '\n\n')
    f.write(str(vars(pr_)) + '\n\n')
    
    f.write('Results:' + '\n\n')
    f.write(str(vars(res_)) + '\n\n')
    
    f.close()
    
    f = open(res_.pic_file, 'w')
    
    pickle.dump(pr_, f)    
    pickle.dump(res_, f)
    
    f.close()
    
def load_results():
    global res_
    global pr_ 
    
    f = open(res_.pic_file, 'r')
    
    pr_ = pickle.load(f)
    res_ = pickle.load(f)
    
    f.close()

    #return prn_, resn_
    
def read_args(args):
    
    for a,b in args.items():
        if not vars(pr_).has_key(a):
            print "Parameter given misspelled"
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
