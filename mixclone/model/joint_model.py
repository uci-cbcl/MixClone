'''
Created on 2013-07-31

@author: Yi Li

pyloh.model.poisson_model

================================================================================

Modified on 2014-04-21

@author: Yi Li
'''
import sys

import numpy as np

from mixclone import constants
from mixclone.preprocess.data import Data
from mixclone.model.model_base import *
from mixclone.model.utils import *

class JointProbabilisticModel(ProbabilisticModel):
    def __init__(self, max_copynumber, subclone_num, baseline_thred):
        ProbabilisticModel.__init__(self, max_copynumber, subclone_num, baseline_thred)
        
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
        self.model_trainer_class = JointModelTrainer


class JointModelTrainer(ModelTrainer):
    def __init__(self, priors, data, max_copynumber, subclone_num, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, max_copynumber, subclone_num, max_iters, stop_value)

    def _init_components(self):
        self.config_parameters = JointConfigParameters(self.max_copynumber, self.subclone_num)
        
        self.model_parameters = JointModelParameters(self.priors, self.data, self.config_parameters)
        
        self.latent_variables = JointLatentVariables(self.data, self.config_parameters)
        
        self.model_likelihood = JointModelLikelihood(self.priors, self.data, self.config_parameters)
    
    #TODO    
    def train(self):
        
        return None
    
    
    def _E_step(self):
        J = self.data.seg_num
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        eta = constants.ETA
        
        rho = self.model_parameters.parameters['rho']
        pi = self.model_parameters.parameters['pi']
        
        for j in range(0, J):
            complete_ll_j = self.model_likelihood.complete_ll_by_seg(self.model_parameters, j)
            complete_ll_j = np.log(rho[j].reshape((1, H))) + np.log(pi.reshape((K, 1))) + complete_ll_j
            
            for h in range(0, H):
                h_T = self.config_parameters.allele_config[h]
            
                if self.data.segments[j].LOH_status == 'FALSE' and check_balance_allele_type(h_T) == False:
                    complete_ll_j[:, h] = np.log(1.0/constants.INF)
            
            psi_j_temp = np.logaddexp.reduce(complete_ll_j, axis=0)
            kappa_j_temp = np.logaddexp.reduce(complete_ll_j, axis=1)
            psi_j = log_space_normalise_rows_annealing(phi_j_temp.reshape((1, H)), eta)
            kappa_j = log_space_normalise_rows_annealing(kappa_j_temp.reshape((1, K)), eta)
            
            self.latent_variables.sufficient_statistics['psi'][j] = psi_j
            self.latent_variables.sufficient_statistics['kappa'][j] = kappa_j
    #TODO        
    def _M_step(self):
        
        return None
    
    def _update_rho(self):
        x = constants.UPDATE_WEIGHTS['x']
        y = constants.UPDATE_WEIGHTS['y']
        
        psi = np.array(self.latent_variables.sufficient_statistics['psi'])
        omega = np.array(self.priors['omega'])
        
        rho_data = psi
        rho_prior = omega/omega.sum()
        
        #rho = x*rho_data + y*rho_priors
        rho = rho_data
        
        self.model_parameters.parameters['rho'] =  rho
    
    def _update_pi(self):
        J = self.data.seg_num
        
        kappa = np.array(self.latent_variables.sufficient_statistics['kappa'])
        
        pi = kappa.sum(axis=0)/J
        
        self.model_parameters.parameters['pi'] = pi


class JointConfigParameters(ConfigParameters):
    def __init__(self, max_copynumber, subclone_num):
        ConfigParameters.__init__(self, max_copynumber, subclone_num)
        
    def _init_components(self):
        self.copynumber = get_copynumber(self.max_copynumber)
        self.copynumber_num = get_copynumber_num(self.max_copynumber)
        self.genotype = get_genotype(self.max_copynumber)
        self.genotype_num = get_genotype_num(self.max_copynumber)
        self.allele_config = get_allele_config(self.max_copynumber)
        self.allele_config_num = get_allele_config_num(self.max_copynumber)
        self.allele_config_CN = get_allele_config_CN(self.max_copynumber)
        self.MU_G = get_MU_G(self.max_copynumber)
        self.Q_HG = get_Q_HG(self.max_copynumber)
        
        
class JointModelParameters(ModelParameters):
    def __init__(self, priors, data, config_parameters):
        ModelParameters.__init__(self, priors, data, config_parameters)
        
    def _init_parameters(self):
        J = self.data.seg_num
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        
        parameters = {}
        parameters['rho'] = np.ones((J, H))*1.0/H
        parameters['pi'] = np.ones(K)*1.0/K
        parameters['phi'] = np.random.random(K)

        self.parameters = parameters
        
    def reinit_parameters(self, phi):
        J = self.data.seg_num
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        
        parameters = {}
        parameters['rho'] = np.ones((J, H))*1.0/H
        parameters['pi'] = np.ones(K)*1.0/K
        parameters['phi'] = phi

        self.parameters = parameters


class JointLatentVariables(LatentVariables):
    def __init__(self, data, config_parameters):
        LatentVariables.__init__(self, data, config_parameters)
        
        self._init_components()
    
    def _init_components(self):
        J = self.data.seg_num
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        
        sufficient_statistics = {}
        sufficient_statistics['psi'] = np.zeros((J, H))
        sufficient_statistics['kappa'] = np.zeros((J, K))
        
        for j in range(0, J):
            sufficient_statistics['psi'][j] = rand_probs(H)
            sufficient_statistics['kappa'][j] = rand_probs(K)
    
        self.sufficient_statistics = sufficient_statistics


class JointModelLikelihood(ModelLikelihood):
    def __init__(self, priors, data, config_parameters):
        ModelLikelihood.__init__(self, priors, data, config_parameters)
        
    def complete_ll_by_seg(self, model_parameters, j):
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        
        ll = np.zeros((K, H))
        
        ll += self._complete_ll_CNA_by_seg(model_parameters, j)
        ll += self._complete_ll_LOH_by_seg(model_parameters, j)
        
        return ll
        
    def _complete_ll_CNA_by_seg(self, model_parameters, j):
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        
        phi = np.array(model_parameters.parameters['phi'])
        
        c_N = constants.COPY_NUMBER_NORMAL
        c_S = constants.COPY_NUMBER_BASELINE
        c_H = np.array(self.config_parameters.allele_config_CN)
        D_N_j = self.data.segments[j].normal_reads_num
        D_T_j = self.data.segments[j].tumor_reads_num
        Lambda_S = self.data.Lambda_S
        
        c_E_j = get_c_E(c_N, c_H, phi)
        lambda_E_j = D_N_j*c_E_j*Lambda_S/c_S
        
        ll_CNA_j = log_poisson_likelihood(D_T_j, lambda_E_j)
        
        return ll_CNA_j
        
    def _complete_ll_LOH_by_seg(self, model_parameters, j):
        H = self.config_parameters.allele_config_num
        K = self.config_parameters.subclone_num
        G = self.config_parameters.genotype_num
        Q_HG = np.array(self.config_parameters.Q_HG).reshape(1, 1, H, G)
        
        phi = np.array(model_parameters.parameters['phi'])
        
        c_N = constants.COPY_NUMBER_NORMAL
        c_H = np.array(self.config_parameters.allele_config_CN)
        mu_N = constants.MU_N
        mu_G = np.array(self.config_parameters.MU_G)
        mu_E = get_mu_E_joint(mu_N, mu_G, c_N, c_H, phi)
        a_T_j = self.data.segments[j].paired_counts[:, 2]
        b_T_j = self.data.segments[j].paired_counts[:, 3]
        d_T_j = a_T_j + b_T_j
        
        ll = np.log(Q_HG) + log_binomial_likelihood_joint(b_T_j, d_T_j, mu_E)
        ll_LOH_j = np.logaddexp.reduce(ll, axis=3).sum(axis=0)
        
        return ll_LOH_j
        
        
        
        
        
