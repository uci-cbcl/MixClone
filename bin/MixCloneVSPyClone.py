#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

def main():
    inbed = open(sys.argv[1])
    inmixclone = open(sys.argv[2])
    inpyclone = open(sys.argv[3])
    
    start_true_lst = []
    end_true_lst = []
    subclone_true_lst = []
    
    start_base_lst = []
    end_base_lst = []
    subclone_base_lst = []
    
    for line in inbed:
        if line[0] == '#':
            continue
        
        fields = line.strip('\n').split('\t')
        
        p_copy = int(fields[3])
        m_copy = int(fields[4])
        
        if p_copy == 1 and m_copy == 1:
            start_base_lst.append(int(fields[1]))
            end_base_lst.append(int(fields[2]))
            subclone_base_lst.append(float(fields[5]))
        else:
            start_true_lst.append(int(fields[1]))
            end_true_lst.append(int(fields[2]))
            subclone_true_lst.append(float(fields[5]))
    
    pos_lst = []
    subclone_pyclone_lst = []
    
    for line in inpyclone:
        if line[0:3] != 'chr':
            continue
        
        fields = line.strip('\n').split('\t')
        pos_lst.append(int(fields[0].split(':')[1]))
        subclone_pyclone_lst.append(float(fields[2]))
        
    start_mixclone_lst = []
    end_mixclone_lst = []
    subclone_mixclone_lst = []
    
    for line in inmixclone:
        if line[0] == '#':
            continue
        
        fields = line.strip('\n').split('\t')
        
        allele_type = fields[10]
        
        if allele_type == 'PM':
            continue
        
        start_mixclone_lst.append(int(fields[2]))
        end_mixclone_lst.append(int(fields[3]))
        subclone_mixclone_lst.append(float(fields[11]))
    
    plt.rcParams['font.size'] = 15
    plt.figure(figsize=(10,10), dpi=150)
    plt.xlim(1, 247249718)
    plt.ylim(0, 1)
    plt.xticks([0, 247249719], [0, 247249719])
    plt.yticks(sp.linspace(0, 1, 11), ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    plt.xlabel('Coordinates of Chromosome 1')
    plt.ylabel('Subclonal cellular prevalence')
    for i in range(0, len(subclone_true_lst)):
        if i == 0:
            plt.plot([start_base_lst[i], end_base_lst[i]], [subclone_base_lst[i]+0.002, subclone_base_lst[i]+0.002], '-', color = 'c', linewidth=1.5, label = 'Ground truth diploid')
        else:
            plt.plot([start_base_lst[i], end_base_lst[i]], [subclone_base_lst[i]+0.002, subclone_base_lst[i]+0.002], '-', color = 'c', linewidth=1.5)
    for i in range(0, len(subclone_true_lst)):
        if i == 0:
            plt.plot([start_true_lst[i], end_true_lst[i]], [subclone_true_lst[i]+0.002, subclone_true_lst[i]+0.002], '-', color = 'r', linewidth=1.5, label = 'Ground truth non-diploid')
        else:
            plt.plot([start_true_lst[i], end_true_lst[i]], [subclone_true_lst[i]+0.002, subclone_true_lst[i]+0.002], '-', color = 'r', linewidth=1.5)
    plt.plot(pos_lst, subclone_pyclone_lst, 'y.', ms=6, label = 'PyClone')
    for i in range(0, len(subclone_mixclone_lst)):
        if i == 0:
            plt.plot([start_true_lst[i], end_true_lst[i]], [subclone_mixclone_lst[i]-0.002, subclone_mixclone_lst[i]-0.002], '-', color = 'b', linewidth=1.5, label = 'MixClone')
        else:
            plt.plot([start_true_lst[i], end_true_lst[i]], [subclone_mixclone_lst[i]-0.002, subclone_mixclone_lst[i]-0.002], '-', color = 'b', linewidth=1.5)
            

    plt.legend(loc=4,prop={'size':10})

    plt.savefig('MixCloneVSPyClone.eps', bbox_inches='tight')
    
    
        
if __name__ == '__main__':
    main()