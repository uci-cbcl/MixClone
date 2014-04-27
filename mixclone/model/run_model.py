'''
Created on 2012-08-13

@author: Yi Li

pyloh.model.run_model

================================================================================

Modified on 2014-04-18

@author: Yi Li
'''
import sys
import time

import numpy as np

from mixclone.model.indep_model import IndepProbabilisticModel
from mixclone.model.joint_model import JointProbabilisticModel

def run_model(args):
    if args.max_copynumber < 1 or args.max_copynumber > 8:
        print "Currently only support max copynumber in range [1, 8]."
        sys.exit(1)
    if args.subclone_num < 1 or args.subclone_num > 5:
        print "Currently only support subclone number in range [1, 5]."
        sys.exit(1)
    
    time_start = time.time()
    
    ## run the independent model
    #indep_model = IndepProbabilisticModel(args.max_copynumber, args.subclone_num, args.baseline_thred)
    #indep_model.read_priors(args.priors)
    #indep_model.read_data(args.filename_base)
    #indep_model.preprocess()
    #indep_model.run(args.max_iters, args.stop_value)
    #indep_model.write_results(args.filename_base)
    
    # run the joint model
    joint_model = JointProbabilisticModel(args.max_copynumber, args.subclone_num, args.baseline_thred)
    joint_model.read_priors(args.priors)
    joint_model.read_data(args.input_filename_base)
    joint_model.preprocess()
    joint_model.run(args.max_iters, args.stop_value)
    joint_model.write_results(args.output_filename_base)
    
    time_end = time.time()
    
    print "*" * 100
    print "* Finish."
    print "* Run time : {0:.2f} seconds".format(time_end - time_start)
    print "* Optimum Log-likelihood : ", joint_model.ll
    print "*" * 100
    sys.stdout.flush()
    
    
    
    