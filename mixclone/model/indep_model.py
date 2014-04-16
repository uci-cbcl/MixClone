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
    
#TODO
class IndepModelTrainer(ModelTrainer):
    def __init__(self, priors, data, max_copynumber, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, max_copynumber, max_iters, stop_value)

    def _init_components(self):
        self.config_parameters = IndepConfigParameters(self.max_copynumber)
        
        self.model_parameters = IndepModelParameters(self.priors, self.data, self.config_parameters)
        
        self.latent_variables = IndepLatentVariables(self.data, self.config_parameters)
        
        self.model_likelihood = IndepModelLikelihood(self.data, self.config_parameters)
        
#TODO    
class IndepConfigParameters(ConfigParameters):
    def __init__(self, max_copynumber):
        ConfigParameters.__init__(self, max_copynumber)
        
    def _init_components(self):
        self.copynumber = get_copynumber(max_copynumber)
        self.copynumber_num = get_copynumber_num(max_copynumber)
        self.allele_config = get_allele_config(max_copynumber)
        self.allele_config_num = get_allele_config_num(max_copynumber)


class IndepModelParameters(ModelParameters):
    def __init__(self, priors, data, config_parameters):
        ModelParameters.__init__(self, priors, data, config_parameters)
        
    def _init_parameters(self):
        J = self.data.seg_num
        
        parameters = {}
        
        parameters['phi'] = np.random.random(J)
        
#TODO        
class IndepLatentVariables(LatentVariables):
    def __init__(self, data, config_parameters):
        LatentVariables.__init__(self, data, config_parameters)
        
        self._init_components()
    
    def _init_components(self):
        
    
    
    