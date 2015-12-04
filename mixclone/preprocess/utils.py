'''
Created on 2013-07-27

@author: Yi Li

pyloh.preprocess.utils

================================================================================

Modified on 2014-04-09

@author: Yi Li
'''
import sys
import numpy as np
from scipy.stats import binom

from mixclone import constants

def BEDParser(bed_file_name):
    inbed = open(bed_file_name)
    
    chroms = []
    starts = []
    ends = []
    
    for line in inbed:
        fields = line.split('\t')
        chrom_name = fields[0]
        chrom_idx = chrom_name_to_idx(chrom_name)
        
        if chrom_idx == -1:
            continue
        
        chrom_name, start, end = fields[0:3]
        
        chroms.append(chrom_name)
        starts.append(int(start))
        ends.append(int(end))
            
    inbed.close()
    
    return (chroms, starts, ends)
    
def chrom_idx_to_name(idx, format):
    if format == 'UCSC':
        chrom_name = 'chr' + str(idx)
    elif format == 'ENSEMBL':
        chrom_name = str(idx)
    else:
        print 'Error: %s not supported' % (format)
        sys.exit(1)

    return chrom_name

def chrom_name_to_idx(chrom_name):
    idx = -1
    
    try:
        idx = int(chrom_name.strip('chr'))
    except:
        pass
        
    return idx

def get_chrom_format(chroms):
    format = 'NONE'
    
    for chrom in chroms:
        if chrom[0:3] == 'chr':
            format = 'UCSC'
            break
        else:
            try:
                idx = int(chrom)
                format = 'ENSEMBL'
                break
            except:
                pass
    
    if format == 'NONE':
        print 'Error: %s not supported' % (chrom)
        sys.exit(-1)
    else:
        return format

def get_chrom_lens_idxs(chrom_idx_list, sam_SQ):
    chrom_lens = []
    chrom_idxs = []
    for i in range(0, len(chrom_idx_list)):
        chrom_idx = chrom_idx_list[i]
    
        for j in range(0, len(sam_SQ)):
            if chrom_idx == chrom_name_to_idx(sam_SQ[j]['SN']):
                chrom_lens.append(int(sam_SQ[j]['LN']))
                chrom_idxs.append(chrom_idx)
                break
        
    return (chrom_lens, chrom_idxs)

def get_segment_name(chrom_name, start, end):
    
    return '_'.join([chrom_name, 'start', str(start), 'end', str(end)])

def normal_heterozygous_filter(counts):
    BAF_N_MAX = constants.BAF_N_MAX
    BAF_N_MIN = constants.BAF_N_MIN
    
    I = counts.shape[0]
    idx_keep = []
    
    for i in xrange(0, I):
        a_N = counts[i, 0]*1.0
        b_N = counts[i, 1]*1.0
        d_N = a_N + b_N
        BAF_N = b_N/d_N
        
        if BAF_N >= BAF_N_MIN and BAF_N <= BAF_N_MAX:
            idx_keep.append(i)
            
    counts = counts[idx_keep]
    
    return counts

def get_BAF_counts(counts):
    BAF_bins = constants.BAF_BINS
    
    a_N = counts[:, 0]*1.0
    b_N = counts[:, 1]*1.0
    a_T = counts[:, 2]*1.0
    b_T = counts[:, 3]*1.0
    
    BAF_N = b_N/(a_N + b_N)
    BAF_T = b_T/(a_T + b_T)
    
    BAF_counts, _, _ = np.histogram2d(BAF_N, BAF_T, bins=(BAF_bins, BAF_bins))
    
    return BAF_counts

def get_LOH_frac(counts):
    I = counts.shape[0]
    
    sites_num_min = constants.SITES_NUM_MIN
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED
  
    if I < sites_num_min:
        LOH_frac = -1
        
        return LOH_frac
        
    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis = 1)
    p_T = binom.cdf(l_T, d_T, p)
    
    LOH_num = np.where(p_T < thred)[0].shape[0]
    LOH_frac = LOH_num*1.0/I
    
    return LOH_frac

def get_LOH_status(LOH_frac, baseline_thred):
    LOH_FRAC_MAX = constants.LOH_FRAC_MAX
    
    if LOH_frac == -1:
        LOH_status = 'NONE'
    elif LOH_frac < baseline_thred:
        LOH_status = 'FALSE'
    elif LOH_frac >= baseline_thred and LOH_frac < LOH_FRAC_MAX:
        LOH_status = 'UNCERTAIN'
    elif LOH_frac >= LOH_FRAC_MAX:
        LOH_status = 'TRUE'
    else:
        LOH_status = 'ERROR'
        
    return LOH_status

def remove_outliers(X):
    std_thred = 0.05
    
    idx_keep = []
    
    n = X.shape[0]
    
    for i in range(0, n):
        if np.abs(X[i] - X.mean()) <= X.std():
            idx_keep.append(i)
    
    if len(idx_keep) == 0 or len(idx_keep) == n:
        return X
    
    X = X[idx_keep]
    
    if X.std() < std_thred:
        return X
    else:
        return remove_outliers(X)


    
