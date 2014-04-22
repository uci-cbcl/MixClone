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

#JointSNVMix
def log_space_normalise_rows_annealing(log_X, eta):
    nrows = log_X.shape[0]
    shape = ( nrows, 1 )
    log_X = log_X*eta
    
    log_norm_const = np.logaddexp.reduce( log_X, axis=1 )
    log_norm_const = log_norm_const.reshape( shape )

    log_X = log_X - log_norm_const
    
    X = np.exp( log_X )
    
    dt = X.dtype
    eps = np.finfo( dt ).eps
    
    X[X <= eps] = 0.
    
    return X

#JointSNVMix
def log_binomial_likelihood(k, n, mu):
    column_shape = (k.size, 1)
    k = k.reshape(column_shape)
    n = n.reshape(column_shape)
    
    row_shape = (1, mu.size)
    mu = mu.reshape(row_shape)
    
    return k * np.log(mu) + (n - k) * np.log(1 - mu)

#JointSNVMix    
def log_binomial_likelihood_multi(k, n, mu):
    d = len(mu.shape)
    column_shape = [k.size]
    column_shape.extend([1 for i in range(0, d)])
    k = k.reshape(column_shape)
    n = n.reshape(column_shape)
    row_shape = [1]
    row_shape.extend(mu.shape)
    mu = mu.reshape(row_shape)
    
    return k * np.log(mu) + (n - k) * np.log(1 - mu)

def log_poisson_likelihood(k, Lambda):
    
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)

def get_omega(max_copynumber):
    tau = constants.TAU
    
    omega = tau*np.ones(max_copynumber+1)
    
    return omega

def get_copynumber(max_copynumber):
    
    return range(0, max_copynumber + 1)
    
def get_copynumber_num(max_copynumber):
    
    return max_copynumber + 1

def get_genotype(max_copynumber):
    genotype = []
    
    for cn in range(0, max_copynumber+1):
        for B_num in range(0, cn+1):
            A_num = cn - B_num
            if A_num == 0 and B_num == 0:
                g_T = 'NULL'
            else:
                g_T = 'A'*A_num + 'B'*B_num
            
            genotype.append(g_T)
    
    return genotype

def get_genotype_num(max_copynumber):
    
    return len(get_genotype(max_copynumber))
    
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
    
def check_balance_allele_type(h_T):
    
    return h_T.count('/') == 0
        
def get_Q_GH(max_copynumber):
    sigma = constants.SIGMA
    
    G = get_genotype_num(max_copynumber)
    H = get_allele_config_num(max_copynumber)
    g_T = get_genotype(max_copynumber)
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
    
def get_c_E(c_N, c_T, phi):
    column_shape = (phi.size, 1)
    row_shape = (1, c_T.size)
    c_T = c_T.reshape(row_shape)
    phi = phi.reshape(column_shape)
    
    return (1 - phi)*c_N + phi*c_T
    
def get_mu_E(mu_N, mu_G, c_N, c_H, phi):
    
    return ((1 - phi)*c_N*mu_N + phi*c_H*mu_G)/((1 - phi)*c_N + phi*c_H)
    
def get_mu_E_multi(mu_N, mu_G, c_N, c_H, phi):
    axis_1_shape = (mu_G.size, 1, 1)
    axis_2_shape = (1, phi.size, 1)
    axis_3_shape = (1, 1, c_H.size)
    mu_G = mu_G.reshape(axis_1_shape)
    phi = phi.reshape(axis_2_shape)
    c_H = c_H.reshape(axis_3_shape)
    
    return ((1 - phi)*c_N*mu_N + phi*c_H*mu_G)/((1 - phi)*c_N + phi*c_H)
    
def rand_probs(N):
    rand_int = np.random.randint(1, N*10, size=N)
    probs = rand_int*1.0/rand_int.sum()
    
    return probs
    
    

