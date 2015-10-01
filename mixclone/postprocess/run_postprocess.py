'''
Created on 2012-08-18

@author: Yi Li

pyloh.postprocess.plot

================================================================================

Modified on 2014-04-11

@author: Yi Li

'''
import os
import sys
import pickle as pkl

import numpy as np
import scipy as sp

PLT_AVAIL = True
try:
    from matplotlib import pyplot as plt
except:
    PLT_AVAIL = False
    
from mixclone import constants
from mixclone.preprocess.data import Data

def run_postprocess(args):
    file_name = args.output_filename_base + '.MixClone.output.pkl'    
    infile = open(file_name, 'rb')
    
    trainer = pkl.load(infile)
    data = trainer.data
    config_parameters = trainer.config_parameters
    model_parameters = trainer.model_parameters
    latent_variables = trainer.latent_variables
    ll = trainer.ll
    
    extract_paired_counts(data, args.output_filename_base)
    extract_segments(data, args.output_filename_base)
    
    if PLT_AVAIL == True:
        extract_BAFheatmap(data, args.output_filename_base)
        extract_seg_plot(data, args.output_filename_base)
    else:
        print "matplotlib.pyplot not available, skip plotting..."
        sys.stdout.flush()
    
    extract_summary(model_parameters, config_parameters,
                    ll, args.output_filename_base)
    
    infile.close()
    
        
def extract_paired_counts(data, output_filename_base):
    counts_file_name = output_filename_base + '.MixClone.counts'
    outfile = open(counts_file_name, 'w')
    segments = data.segments
    
    print "Extracting paired counts file..."
    sys.stdout.flush()
    
    outfile.write('\t'.join(['#seg_name', 'normal_A', 'normal_B', 'tumor_A',
                             'tumor_B', 'chrom', 'pos']) + '\n')
    
    for j in range(0, data.seg_num):
        for i in range(0, segments[j].paired_counts.shape[0]):
            outfile.write(segments[j].name + '\t'
                          + '\t'.join(map(str, segments[j].paired_counts[i])) + '\n')

    outfile.close()
    
    
def extract_segments(data, output_filename_base):
    segments_file_name = output_filename_base + '.MixClone.segments'
    outfile = open(segments_file_name, 'w')
    segments = data.segments
    
    print "Extracting segments file..."
    sys.stdout.flush()
    
    outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num',
                             'tumor_reads_num', 'LOH_frac', 'LOH_status', 'log2_ratio',
                             'copy_number', 'allele_type', 'subclone_prev', 'subclone_cluster']) + '\n')
    
    for j in range(0, data.seg_num):
        outfile.write('\t'.join(map(str, [segments[j].name, segments[j].chrom_name,segments[j].start,
                                segments[j].end, segments[j].normal_reads_num, segments[j].tumor_reads_num,
                                segments[j].LOH_frac, segments[j].LOH_status, segments[j].log2_ratio,
                                segments[j].copy_number, segments[j].allele_type,
                                "{0:.3f}".format(segments[j].subclone_prev), segments[j].subclone_cluster])) + '\n')

    outfile.close()
    

def extract_BAFheatmap(data, output_filename_base):
    BAF_counts_min = constants.BAF_COUNTS_MIN
    BAF_counts_max = constants.BAF_COUNTS_MAX    

    outheatmap_dir_name = output_filename_base + '.MixClone.heatmap'
    if os.path.exists(outheatmap_dir_name) == False:
        os.mkdir(outheatmap_dir_name)
        
    seg_num = data.seg_num
    
    for j in range(0, seg_num):
        BAF_counts_j = data.segments[j].BAF_counts
        seg_name_j = data.segments[j].name
        BAF_counts_sub = BAF_counts_j[BAF_counts_min:BAF_counts_max, BAF_counts_min:BAF_counts_max]
        color_max_j = BAF_counts_sub.max()
        
        print 'Plotting segment {0}...'.format(seg_name_j)
        sys.stdout.flush()
        
        plt.figure(figsize=(8,8), dpi = 150)
        plt.xlim((0, 100))
        plt.ylim((0, 100))
        plt.xticks(sp.linspace(0, 100, 11), sp.linspace(0, 1, 11))
        plt.yticks(sp.linspace(0, 100, 11), sp.linspace(0, 1, 11))
        plt.xlabel('Tumor sample B allele frequency')
        plt.ylabel('Normal sample B allele frequency')
        plt.imshow(BAF_counts_j, vmin = 0, vmax = max(1, color_max_j))
        cbar = plt.colorbar(ticks=[0, color_max_j], orientation='vertical', shrink=0.78)
        cbar.ax.set_yticklabels(['0', '>= ' + str(int(color_max_j))])
        plt.savefig('./' + outheatmap_dir_name + '/' + seg_name_j, bbox_inches='tight')


def extract_seg_plot(data, output_filename_base):
    subclone_plot_file_name = output_filename_base + '.MixClone.segplot.png'
    
    subclone_prev_lst = []
    copynumber_lst = []
    seg_num = data.seg_num

    print "Extracting segments plot file..."
    sys.stdout.flush()

    for j in range(0, seg_num):
        if data.segments[j].baseline_label == True or data.segments[j].allele_type == 'PM':
            continue
        subclone_prev_lst.append(data.segments[j].subclone_prev)
        copynumber_lst.append(data.segments[j].copy_number)

    X = len(subclone_prev_lst)
    
    plt.figure(figsize=(8,8), dpi = 150)
    plt.plot(range(1, X+1), subclone_prev_lst, 'o')
    plt.xlim(0, X+1)
    plt.ylim(0, 1)
    plt.xlabel('Copy number')
    plt.ylabel('Subclonal cellular prevalence')
    plt.xticks(sp.linspace(1, X, X), copynumber_lst)
    plt.yticks(sp.linspace(0, 1, 11), ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    plt.savefig(subclone_plot_file_name, bbox_inches='tight')

    
def extract_summary(model_parameters, config_parameters, ll, output_filename_base):
    summary_file_name = output_filename_base + '.MixClone.summary'
    outfile = open(summary_file_name, 'w')
    
    phi = model_parameters.parameters['phi']
    subclone_cluster = '\t'.join(map(str, range(1, config_parameters.subclone_num+1)))
    subclone_prev = '\t'.join(map("{0:.3f}".format, phi.tolist()))

    print "Extracting summary file..."
    sys.stdout.flush()
    
    outfile.write("Model : joint\n")
    outfile.write("Maximum copy number : %s\n" % (config_parameters.max_copynumber))
    outfile.write("Subclone number : %s\n" % (config_parameters.subclone_num))
    outfile.write("Subclone cluster :             %s\n" % subclone_cluster)
    outfile.write("Subclonal cellular prevalence : %s\n" % subclone_prev)
    outfile.write("Optimum log-likelihood : %s\n" % (ll))
    
    outfile.close()
    
