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

def get_copynumber(max_copynumber):
    
    return range(0, max_copynumber + 1)
    
def get_copynumber_num(max_copynumber):
    
    return max_copynumber + 1

def get_genotypes(max_copynumber):
    genotypes = []
    
    for cn in range(0, max_copynumber+1):
        for B_num in range(0, cn+1):
            A_num = cn - B_num
            if A_num == 0 and B_num == 0:
                g_T = 'NULL'
            else:
                g_T = 'A'*A_num + 'B'*B_num
            
            genotypes.append(g_T)
    
    return genotypes

def get_genotypes_num(max_copynumber):
    
    return len(get_genotypes(max_copynumber))
    
def get_MU_G(max_copynumber):
    empiri_BAF = constants.EMPIRI_BAF
    empiri_AAF = constants.EMPIRI_AAF
    err = constants.ERR
    
    MU_G = []
    
    for cn in range(0, max_copynumber+1):
        for B_num in range(0, cn+1):
            A_num = cn - B_num
            if A_num == 0 and B_num == 0:
                mu_G = empiri_BAF/(empiri_BAF + empiri_AAF)
            elif A_num == 0 and B_num != 0:
                mu_G = 1 - err
            elif A_num != 0 and B_num == 0:
                mu_G = err
            else:
                mu_G = empiri_BAF*B_num/(empiri_BAF*B_num + empiri_AAF*A_num)
                
            MU_G.append(mu_G)
    
    return MU_G



