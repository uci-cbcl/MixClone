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

def get_allele_config(max_copynumber):
    allele_config = []
    
    for cn in range(0, max_copynumber+1):
        for M_num in range(0, (cn+2)/2):
            P_num = cn - M_num
                        
            if P_num == 0 and M_num == 0:
                h_T = 'NULL'
            elif P_num == M_num:
                h_T = 'P'*P_num + 'M'*M_num
            else:
                h_T = 'P'*P_num + 'M'*M_num + '/' + 'P'*M_num + 'M'*P_num
                
            allele_config.append(h_T)
            
    return allele_config
    
def get_allele_config_num(max_copynumber):
    
    return len(get_allele_config(max_copynumber))
    
def get_allele_config_CN(max_copynumber):
    allele_config_CN = []
    
    for cn in range(0, max_copynumber+1):
        for M_num in range(0, (cn+2)/2):
            
            allele_config_CN.append(cn)
            
    return allele_config_CN

def check_GH_compat(g, h):
    if g == 'NULL' and 'h' == 'NULL':
        return True
    
    h_half = h.split('/')[0]
    
    P_num = h_half.count('P')
    M_num = h_half.count('M')
    A_num = g.count('A')
    B_num = g.count('B')
    
    if A_num == P_num and B_num == M_num:
        return True
    elif A_num == M_num and B_num == P_num:
        return True
    else:
        return False
        
def get_Q_GH(max_copynumber):
    sigma = constants.SIGMA
    
    G = get_genotypes_num(max_copynumber)
    H = get_allele_config_num(max_copynumber)
    g_T = get_genotypes(max_copynumber)
    h_T = get_allele_config(max_copynumber)
    
    Q_GH = np.ones((G, H))*sigma
    
    for h in range(0, H):
        if h_T[h].count('/') == 0:
            compat_num = 1
        elif h_T[h].count('/') == 1:
            compat_num = 2
        
        for g in range(0, G):
            if check_GH_compat(g_T[g], h_T[h]) == True:
                Q_GH[g, h] = (1 - sigma*(G - compat_num ))/compat_num 
            
    return Q_GH

