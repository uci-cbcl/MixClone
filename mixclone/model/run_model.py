'''
Created on 2012-08-13

@author: Yi Li

pyloh.model.run_model

================================================================================

Modified on 2014-04-18

@author: Yi Li
'''
import sys

import numpy as np

from mixclone.model.indep_model import IndepProbabilisticModel

def run_model(args):
    # run the independent model
    indep_model = IndepProbabilisticModel(args.max_copynumber, args.subclone_num, args.baseline_thred)
    indep_model.read_priors(args.priors)
    indep_model.read_data(args.filename_base)
    indep_model.preprocess()
    indep_model.run(args.max_iters, args.stop_value)
    indep_model.write_results(args.filename_base)
    
    print "*" * 100
    print "* Finish."
    print "*" * 100
    sys.stdout.flush()
    
    
    
    