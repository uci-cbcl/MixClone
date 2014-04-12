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

parser_preprocess.add_argument('filename_base',
                          help='''Base name of preprocessed files to be created.''')

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
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)
