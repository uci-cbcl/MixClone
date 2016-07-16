'''
Created on 2013-07-20

@author: Yi Li

pyloh.constants

================================================================================

Modified on 2014-04-05

@author: Yi Li
'''
import numpy as np

BAF_N_MIN = 0.4
BAF_N_MAX = 0.6
BAF_COUNTS_MIN = 10
BAF_COUNTS_MAX = 95
BAF_BINS = np.array(range(0, 100 + 1))/100.0
LOH_FRAC_MAX = 0.25
SITES_NUM_MIN = 20
PHI_INIT = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
PHI_RANGE = [i/100.0 for i in range(2, 99)]
BINOM_TEST_P = 0.5
BINOM_TEST_THRED = 0.025

BURN_IN = 10
EPS = np.finfo(float).eps
ETA = 1.0
INF = float('inf')
LL_RATIO_CHANGE_THRED = 0.01
LL_PERCENT_CHANGE_THRED = 0.90

TAU = 500
SIGMA = 0.0001
ERR = 0.01
EMPIRI_BAF = 0.485
EMPIRI_AAF = 1 - EMPIRI_BAF

MU_N = EMPIRI_BAF/(EMPIRI_BAF + EMPIRI_AAF)
COPY_NUMBER_NORMAL = 2
COPY_NUMBER_BASELINE = 2
ALLELE_TYPE_BASELINE = 'PM'

UPDATE_WEIGHTS = {}
UPDATE_WEIGHTS['x'] = 0.9
UPDATE_WEIGHTS['y'] = 0.1



CHROM_START = 0

CHROM_IDX_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
              20, 21, 22]
