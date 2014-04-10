'''
Created on 2013-08-13

@author: Yi Li

pyloh.preprocess.data

================================================================================

Modified on 2014-04-09

@author: Yi Li
'''
import sys
import pickle
import numpy as np

from mixclone import constants
from mixclone.preprocess.utils import *


class Segment:
    def __init__(self):
        self.name = ""
        self.chrom_idx = -1
        self.starts = -1 
        self.ends = -1
        self.normal_reads_num = -1
        self.tumor_reads_num = -1
        self.sites_num = 0
        self.LOH_frac = 0.0
        self.LOH_status = 'NONE'
        self.log2_ratio = 0.0
        self.paired_counts = []
        self.BAF_counts = None
        self.copy_number = -1
        self.allele_type_idx = ""
        
        
        