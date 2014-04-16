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
        
    