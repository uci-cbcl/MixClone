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
        self.chrom_name = ""
        self.start = -1 
        self.end = -1
        self.normal_reads_num = -1
        self.tumor_reads_num = -1
        self.sites_num = 0
        self.LOH_frac = 0.0
        self.LOH_status = 'NONE'
        self.baseline_label = 'FALSE'
        self.log2_ratio = 0.0
        self.paired_counts = None
        self.BAF_counts = None
        self.copy_number = -1
        self.allele_type = 'NONE'
        self.subclone_prev = -1
        self.subclone_cluster = 'NONE'
        

class Data:
    def __init__(self):
        self.seg_num = 0
        self.Lambda_S = -1
        self.segments = []
        
    def load_segments(self, normal_bam, tumor_bam, bed_file_name):
        chrom_idx_list = constants.CHROM_IDX_LIST
        chrom_start = constants.CHROM_START
        
        sam_SQ = normal_bam.header['SQ']
        sam_chrom_format = get_chrom_format(map(lambda x:x['SN'], sam_SQ))
        chrom_lens, chrom_idxs = get_chrom_lens_idxs(chrom_idx_list, sam_SQ)
        
        bed_chroms, bed_starts, bed_ends = BEDParser(bed_file_name)
        bed_chrom_format = get_chrom_format(bed_chroms)
        bed_num = len(bed_chroms)
        
        for i in range(0, bed_num):
            chrom_idx = chrom_name_to_idx(bed_chroms[i])
            chrom_name = chrom_idx_to_name(chrom_idx, sam_chrom_format)
            seg_name = get_segment_name(chrom_name, bed_starts[i], bed_ends[i])
            
            if chrom_idx not in chrom_idx_list:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue
            
            chrom_lst_idx = chrom_idxs.index(chrom_idx)
            
            if bed_starts[i] < chrom_start or bed_ends[i] > chrom_lens[chrom_lst_idx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue

            normal_reads_num = normal_bam.count(chrom_name, bed_starts[i], bed_ends[i])
            tumor_reads_num = tumor_bam.count(chrom_name, bed_starts[i], bed_ends[i])
            
            segment_i = Segment()
            segment_i.name = seg_name
            segment_i.chrom_idx = chrom_idx
            segment_i.chrom_name = chrom_name
            segment_i.start = bed_starts[i]
            segment_i.end = bed_ends[i]
            segment_i.normal_reads_num = normal_reads_num
            segment_i.tumor_reads_num = tumor_reads_num
            segment_i.log2_ratio = np.log2(1.0*tumor_reads_num/normal_reads_num)
            
            self.segments.append(segment_i)
            self.seg_num += 1
    
    def get_LOH_frac(self):        
        for j in range(0, self.seg_num):
            self.segments[j].LOH_frac = get_LOH_frac(self.segments[j].paired_counts)
            
    def get_LOH_status(self, baseline_thred):
        for j in range(0, self.seg_num):
            self.segments[j].LOH_status = get_LOH_status(self.segments[j].LOH_frac, baseline_thred)
    
    def compute_Lambda_S(self):
        reads_depth_ratio = []
        
        for j in range(0, self.seg_num):
            if self.segments[j].LOH_status == 'FALSE':
                ratio = self.segments[j].tumor_reads_num*1.0/self.segments[j].normal_reads_num
                reads_depth_ratio.append(ratio)
                
        reads_depth_ratio = np.array(reads_depth_ratio)
        reads_depth_ratio = remove_outliers(reads_depth_ratio)
        
        for j in range(0, self.seg_num):
            ratio = self.segments[j].tumor_reads_num*1.0/self.segments[j].normal_reads_num
            
            if self.segments[j].LOH_status == 'FALSE' and ratio in reads_depth_ratio:
                self.segments[j].baseline_label = 'TRUE'
        
        if reads_depth_ratio.shape[0] != 0:
            self.Lambda_S = reads_depth_ratio.mean()
        else:
            print 'Error: No diploid segments found, existing...'
            sys.exit(1)
    

