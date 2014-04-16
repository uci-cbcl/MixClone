'''
Created on 2010-08-30 for log_binomial_likelihood
Created on 2010-11-22 for log_space_normalise_rows_annealing

@author: Andrew Roth

JointSNVMix-0.6.2
joint_snv_mix.classification.utils.log_pdf.log_binomial_likelihood
joint_snv_mix.classification.utils.normalise.log_space_normalise_rows

================================================================================

Modified on 2013-07-31

@author: Yi Li

pyloh.model.utils

================================================================================

Modified on 2014-04-15

@author: Yi Li
'''
import numpy as np
from scipy.special import gammaln

from mixclone import constants

def get_omega(max_copynumber):
    tau = constants.TAU
    
    omega = tau*np.ones(max_copynumber+1)
    
    return omega

def get_copynumber_tumor(max_copynumber):
    
    return range(0, max_copynumber + 1)
    
def get_copynumber_tumor_num(max_copynumber):
    
    return max_copynumber + 1