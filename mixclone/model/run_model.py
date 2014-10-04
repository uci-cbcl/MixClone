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
from mixclone.model.utils import model_selection_by_ll

def run_model(args):
    if args.max_copynumber < 1 or args.max_copynumber > 8:
        print "Currently only support max copynumber in range [1, 8]."
        sys.exit(1)
    if args.subclone_num < 0:
        print "Subclone number must be positive integers."
        sys.exit(1)
        
    if args.subclone_num == 0:
        run_all_subclone(args)
    else:
        run_one_subclone(args)
    
        
def run_one_subclone(args):
    time_start = time.time()
    
    # run the joint model
    joint_model = JointProbabilisticModel(args.max_copynumber, args.subclone_num, \
                                          args.baseline_thred)
    joint_model.read_data(args.input_filename_base)
    joint_model.preprocess()
    joint_model.run(args.max_iters, args.stop_value)
    joint_model.write_results(args.output_filename_base)
    
    time_end = time.time()
    
    print "*" * 100
    print "* Finish."
    print "* Run time : {0:.2f} seconds".format(time_end - time_start)
    print "* Optimum log-likelihood : ", joint_model.trainer.ll
    print "*" * 100
    sys.stdout.flush()
    

def run_all_subclone(args):
    time_start = time.time()
    
    ll_lst = []
    subclone_num_lst = []
    
    for subclone_num in range(1, 6):
        # run the joint model
        output_filename_base_k = args.output_filename_base + '_subclone_num_' + \
                                 str(subclone_num)
        
        joint_model = JointProbabilisticModel(args.max_copynumber, subclone_num,
                                              args.baseline_thred)
        joint_model.read_data(args.input_filename_base)
        joint_model.preprocess()
        joint_model.run(args.max_iters, args.stop_value)
        joint_model.write_results(output_filename_base_k)
        
        ll_lst.append(joint_model.trainer.ll)
        subclone_num_lst.append(subclone_num)
    
    time_end = time.time()
    
    run_time = time_end - time_start
    
    get_summary(ll_lst, subclone_num_lst, run_time, args.output_filename_base)


def get_summary(ll_lst, subclone_num_lst, run_time, output_filename_base):
    model_selection_filename = output_filename_base + '.MixClone.model_selection'
    outfile = open(model_selection_filename, 'w')

    subclone_num_optimum, ll_change_ratio, ll_change_percent, ll_change_ratio_total = \
    model_selection_by_ll(ll_lst, subclone_num_lst)
    
    print "*" * 100
    print "* Finish."
    print "* Run time : {0:.2f} seconds".format(run_time)
    print "*" * 100
    
    for i in range(0, 5):
        print "* Log-likelihood for subclone number %s : %s" % (i+1, ll_lst[i])
        outfile.write('Log-likelihood for subclone number %s : %s\n' % (i+1, ll_lst[i]))
    
    print "*" * 100
    outfile.write("*" * 100 + '\n')
        
    for i in range(0, 4):
        print "* Log-likelihood change ratio for subclone number %s -> %s : %s" \
        % (i+1, i+2, ll_change_ratio[i])
        outfile.write('Log-likelihood change ratio for subclone number %s -> %s : %s\n' \
        % (i+1, i+2, ll_change_ratio[i]))
        
    print "* Log-likelihood change ratio in total for subclone number %s -> %s : %s" \
    % (1, 5, ll_change_ratio_total)
    outfile.write('Log-likelihood change ratio in total for subclone number %s -> %s : %s\n' \
    % (1, 5, ll_change_ratio_total))

    print "*" * 100
    outfile.write("*" * 100 + '\n')

    for i in range(0, 4):
        print "* Log-likelihood cumulative change percent for subclone number %s -> %s : %s" \
        % (1, i+2, ll_change_percent[i])
        outfile.write('Log-likelihood cumulative change percent for subclone number %s -> %s : %s\n' \
        % (1, i+2, ll_change_percent[i]))

    print "*" * 100
    outfile.write("*" * 100 + '\n')
    
    print "* Suggested subclone number : ", subclone_num_optimum
    print "*" * 100
    sys.stdout.flush()
    
    outfile.write('Suggested subclone number : %s\n' % (subclone_num_optimum))
    
    outfile.close()
