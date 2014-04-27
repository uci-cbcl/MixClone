#!/usr/bin/env python

#=======================================================================================================================
# Created on 2014-04-11
# @author: Yi Li
#
# MixClone
# Copyright (c) 2014 Yi Li <yil8@uci.edu>
#
# This code is free software; you can redistribute it and/or modify it
# under the terms of GNU GPL v2.0 (see the file LICENSE included with the distribution).

# Most features are extended from PyLOH-1.2.0 (https://github.com/uci-cbcl/PyLOH)
# Some features are built on top of JointSNVMix-0.6.2 (http://code.google.com/p/joint-snv-mix/).
#=======================================================================================================================
import argparse

from mixclone.preprocess.run_preprocess import run_preprocess
from mixclone.model.run_model import run_model
from mixclone.postprocess.run_postprocess import run_postprocess


parser = argparse.ArgumentParser(prog='MixClone')
subparsers = parser.add_subparsers()

#===============================================================================
# Add preprocess sub-command
#===============================================================================
parser_preprocess = subparsers.add_parser('preprocess',
                                    help='''Preprocess paired normal and tumor BAM files''')

parser_preprocess.add_argument('reference_genome',
                          help='''FASTA file for reference genome.''')

parser_preprocess.add_argument('segments_bed',
                          help='''BED file for segments.''')

parser_preprocess.add_argument('normal_bam',
                          help='''BAM file for normal sample.''')

parser_preprocess.add_argument('tumor_bam',
                          help='''BAM file for tumor sample.''')

parser_preprocess.add_argument('input_filename_base',
                          help='''Base name of the preprocessed input file to be created.''')

parser_preprocess.add_argument('--min_depth', default=20, type=int,
                          help='''Minimum reads depth required for both normal and tumor samples. 
                          Default is 20.''')

parser_preprocess.add_argument('--min_base_qual', default=10, type=int,
                          help='''Minimum base quality required. Default is 10.''')

parser_preprocess.add_argument('--min_map_qual', default=10, type=int,
                          help='''Minimum mapping quality required. Default is 10.''')

parser_preprocess.add_argument('--process_num', default=1, type=int,
                          help='''Number of processes to launch for preprocessing. Default is 1.''')

parser_preprocess.set_defaults(func=run_preprocess)

#===============================================================================
# Add run_model sub-command
#===============================================================================
parser_run_model = subparsers.add_parser('run_model',
                                      help='''Run a probabilistic model based analysis. Requires preprocessed
                                      input file that have been created.''')

parser_run_model.add_argument('input_filename_base',
                            help='Base name of the preprocessed input file created.')

parser_run_model.add_argument('output_filename_base',
                            help='Base name of the output file to be created.')

parser_run_model.add_argument('--max_copynumber', default=4, type=int,
                            help='''Maximum copy number of each segment allows to take. Default is 4.''')

parser_run_model.add_argument('--subclone_num', default=-1, type=int,
                            help='''Number of subclones within the tumor sample. If not provided,
                                go through [1, 5] and select the best model. Default is None.''')

parser_run_model.add_argument('--baseline_thred', default=0.09, type=float,
                            help='''The threshold of LOH SNP sites fraction within each segment to
                            define the segment as baseline. Default is 0.09.''')

parser_run_model.add_argument('--priors', default=None, type=str,
                             help='''File of the prior distribution. If not provided,
                                use uniform prior. Default is None.''')

parser_run_model.add_argument('--max_iters', default=100, type=int,
                          help='''Maximum number of iterations for training. Default is 100.''')

parser_run_model.add_argument('--stop_value', default=1e-7, type=float,
                          help='''Stop value of the EM algorithm for training. If the change of log-likelihood is lower
                          than this value, stop training. Default is 1e-7.''')

parser_run_model.set_defaults(func=run_model)

#===============================================================================
# Add postprocess sub-command
#===============================================================================
parser_postprocess = subparsers.add_parser('postprocess',
                                      help='''Extract various result files from the outputfile.''')

parser_postprocess.add_argument('output_filename_base',
                            help='''Base name of the output file created.''')

parser_postprocess.add_argument('--input', default=False, action='store_true',
                          help='''Extract from *.Mixclone.input.pkl instead of *.Mixclone.output.pkl.
                          Default is False.''')

parser_postprocess.set_defaults(func=run_postprocess)

#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)
