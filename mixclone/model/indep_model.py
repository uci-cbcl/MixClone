'''
Created on 2013-07-31

@author: Yi Li

pyloh.model.poisson_model

================================================================================

Modified on 2014-04-15

@author: Yi Li
'''
import sys

import numpy as np

from mixclone import constants
from mixclone.preprocess.data import Data
from mixclone.model.model_base import *
from mixclone.model.utils import *

class IndepProbabilisticModel(ProbabilisticModel):
    def __init__(self, max_copynumber, baseline_thred):
        ProbabilisticModel.__init__(self, max_copynumber, baseline_thred)
        
    def read_priors(self, priors_filename):
        if priors_filename != None:
            self.priors_parser.read_priors(priors_filename, self.max_copynumber)
            self.priors = self.priors_parser.priors
        else:
            self.priors = {}
            self.priors['omega'] = np.array(get_omega(self.max_copynumber))*1.0
            
    def preprocess(self):
        #self.data.get_LOH_frac()
        self.data.get_LOH_status(self.baseline_thred)
        self.data.compute_Lambda_S()
        
    def _init_components(self):
        self.model_trainer_class = IndepModelTrainer
    

class IndepModelTrainer(ModelTrainer):
    def __init__(self, priors, data, max_copynumber, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, max_copynumber, max_iters, stop_value)

    def _init_components(self):
        self.config_parameters = IndepConfigParameters(self.max_copynumber)
        
        self.model_parameters = IndepModelParameters(self.priors, self.data, self.config_parameters)
        
        self.latent_variables = IndepLatentVariables(self.data, self.config_parameters)
        
        self.model_likelihood = IndepModelLikelihood(self.priors, self.data, self.config_parameters)
    
    def train(self):
        seg_num = self.data.seg_num
        
        for j in range(0, seg_num):
            if self.data.segments[j].LOH_status == 'NONE':
                continue
            elif self.data.segments[j].LOH_status == 'FALSE':#bug
                self.data.segments[j].allele_type = constants.ALLELE_TYPE_BASELINE
                self.data.segments[j].copy_number = constants.COPY_NUMBER_BASELINE
                continue
            
            h_j, c_H_j, phi_j = self.train_by_seg(j)
            
            self.data.segments[j].allele_type = h_j
            self.data.segments[j].copy_number = c_H_j
            self.data.segments[j].subclone_prev = phi_j
            
    def train_by_seg(self, j):
        H = self.config_parameters.allele_config_num
        ll_lst = []
        phi_lst = []
        
        for h in range(0, H):
            ll, phi = self.bisec_search_ll(j, h)
            ll_lst.append(ll)
            phi_lst.append(phi)
            
            h_ = self.config_parameters.allele_config[h]
            c_H_ = self.config_parameters.allele_config_CN[h]
            
            self._print_running_info(j, h_, c_H_, phi, ll)
            
        ll_lst = np.array(ll_lst)
        idx_optimum = ll_lst.argmax()
        
        h_optimum = self.config_parameters.allele_config[idx_optimum]
        c_H_optimum = self.config_parameters.allele_config_CN[idx_optimum]
        phi_optimum = phi_lst[idx_optimum]
        
        return (h_optimum, c_H_optimum, phi_optimum)
        
    def bisec_search_ll(self, j, h):
        phi_start = 0.01
        phi_end = 0.99
        phi_stop = 1e-5
        phi_change = 1
        
        while phi_change > phi_stop:
            phi_left = phi_start + (phi_end - phi_start)*1/3
            phi_right = phi_start + (phi_end - phi_start)*2/3
            
            ll_left = self.model_likelihood.ll_by_seg(h, phi_left, j)
            ll_right = self.model_likelihood.ll_by_seg(h, phi_right, j)
            
            if ll_left >= ll_right:
                phi_change = phi_end - phi_right
                phi_end = phi_right
            else:
                phi_change = phi_left - phi_start
                phi_start = phi_left
                            
        phi_optimum = (phi_start + phi_end)/2
        ll_optimum = self.model_likelihood.ll_by_seg(h, phi_optimum, j)
        
        return (ll_optimum, phi_optimum)

    def _print_running_info(self, j, h_j, c_H_j, phi_j, ll_j):
        print "#" * 100
        print "# Running Info."
        print "#" * 100
        print "Model : independent"
        print "Segment : ", self.data.segments[j].name
        print "Estimated copy number: ", c_H_j
        print "Estimated allele type : ", h_j
        print "Estimated subclone cellular prevalence : ", phi_j
        print "Log-likelihood : ", ll_j
        sys.stdout.flush()
        

class IndepConfigParameters(ConfigParameters):
    def __init__(self, max_copynumber):
        ConfigParameters.__init__(self, max_copynumber)
        
    def _init_components(self):
        self.copynumber = get_copynumber(self.max_copynumber)
        self.copynumber_num = get_copynumber_num(self.max_copynumber)
        self.genotype = get_genotype(self.max_copynumber)
        self.genotype_num = get_genotype_num(self.max_copynumber)
        self.allele_config = get_allele_config(self.max_copynumber)
        self.allele_config_num = get_allele_config_num(self.max_copynumber)
        self.allele_config_CN = get_allele_config_CN(self.max_copynumber)
        self.MU_G = get_MU_G(self.max_copynumber)
        self.Q_GH = get_Q_GH(self.max_copynumber)


class IndepModelParameters(ModelParameters):
    def __init__(self, priors, data, config_parameters):
        ModelParameters.__init__(self, priors, data, config_parameters)
        
    def _init_parameters(self):
        J = self.data.seg_num
        
        parameters = {}
        parameters['phi'] = np.random.random(J)
        
        self.parameters = parameters
        
        
class IndepLatentVariables(LatentVariables):
    def __init__(self, data, config_parameters):
        LatentVariables.__init__(self, data, config_parameters)
        
        self._init_components()
    
    def _init_components(self):
        J = self.data.seg_num
        H = self.config_parameters.allele_config_num
        
        latent_variables = {}
        latent_variables['H'] = np.random.randint(0, H, J)
        
        self.latent_variables = latent_variables


class IndepModelLikelihood(ModelLikelihood):
    def __init__(self, priors, data, config_parameters):
        ModelLikelihood.__init__(self, priors, data, config_parameters)
    
    def ll_by_seg(self, h, phi, j):
        ll = 0
        
        ll += self._ll_CNA_by_seg(h, phi, j)
        ll += self._ll_LOH_by_seg(h, phi, j)
        
        return ll
        
    def _ll_CNA_by_seg(self, h, phi, j):
        c_N = constants.COPY_NUMBER_NORMAL
        c_S = constants.COPY_NUMBER_BASELINE
        c_H = self.config_parameters.allele_config_CN[h]
        D_N_j = self.data.segments[j].normal_reads_num
        D_T_j = self.data.segments[j].tumor_reads_num
        Lambda_S = self.data.Lambda_S
        
        c_E_j = get_c_E(c_N, c_H, phi)
        lambda_E_j = np.array(D_N_j*c_E_j*Lambda_S/c_S)
        
        ll_CNA_j = log_poisson_likelihood(D_T_j, lambda_E_j).sum()
        
        return ll_CNA_j
    
    def _ll_LOH_by_seg(self, h, phi, j):
        Q_GH = np.array(self.config_parameters.Q_GH)
        eta = constants.ETA
        c_N = constants.COPY_NUMBER_NORMAL
        c_H = self.config_parameters.allele_config_CN[h]
        mu_N = constants.MU_N
        mu_G = np.array(self.config_parameters.MU_G)
        mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
        a_T_j = self.data.segments[j].paired_counts[:, 2]
        b_T_j = self.data.segments[j].paired_counts[:, 3]
        d_T_j = a_T_j + b_T_j
        
        ll = np.log(Q_GH[:, h]) + log_binomial_likelihood(b_T_j, d_T_j, mu_E)
        
        ll_LOH_j = np.logaddexp.reduce(ll, axis=1).sum()
        
        return ll_LOH_j
    
    
    