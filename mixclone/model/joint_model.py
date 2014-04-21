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


class JointConfigParameters(ConfigParameters):
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
        
        




