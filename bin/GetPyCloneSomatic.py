#!/usr/bin/env python

import sys
from collections import Counter

import pysam
import numpy as np

def reads_to_bases(reads, pos):
    bases = []
    
    for read in reads:
        if pos not in read.positions:
            continue
        
        idx = read.positions.index(pos)
        base = read.seq[idx-1:idx].upper()
        bases.append(base)
    
    return bases

def get_CN(pos, chrom_vcf, starts, ends, chroms, p_copy, m_copy):
    J = len(starts)
    
    minor_cn = -1
    major_cn = -1

    for j in range(0, J):
        if pos >= starts[j] and pos < ends[j]:
            minor_cn = min(p_copy[j], m_copy[j])
            major_cn = max(p_copy[j], m_copy[j])
            
            return (minor_cn, major_cn)
        
    print '%s\t%s not found in segments...' % (chrom_vcf, pos)
        
        

def main():
    innormal = sys.argv[1]
    intumor = sys.argv[2]
    infasta = sys.argv[3]
    invcf = sys.argv[4]
    inseg = sys.argv[5]
    outcounts = sys.argv[6]
    chrom = sys.argv[7]
    
    normal_bam = pysam.Samfile(innormal, 'rb')
    tumor_bam = pysam.Samfile(intumor, 'rb')
    ref_genome_fasta = pysam.Fastafile(infasta)
    invcf = open(invcf)
    inseg = open(inseg)
    outcounts = open(outcounts, 'w')
    
    normal_pileup = normal_bam.pileup(chrom)
    tumor_pileup = tumor_bam.pileup(chrom)
    
    starts = []
    ends = []
    chroms = []
    p_copy = []
    m_copy = []
    
    for line in inseg:
        if line[0] == '#':
            continue
        
        fields = line.strip('\n').split('\t')
        chroms.append(fields[0])
        starts.append(int(fields[1]))
        ends.append(int(fields[2]))
        p_copy.append(int(fields[3]))
        m_copy.append(int(fields[4]))
        
    outcounts.write('\t'.join(['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']) + '\n')
    for line in invcf:
        if line[0] == '#':
            continue
        
        chrom_vcf, pos, __, ref, alt = line.strip('\n').split('\t')[0:5]
        pos = int(pos)
        ref = ref.upper()
        alt = alt.upper()
        
        if chrom_vcf != chrom:
            continue
        
        reads_N = normal_bam.fetch(chrom, pos - 1, pos)
        reads_T = tumor_bam.fetch(chrom, pos - 1, pos)
        
        bases_N = reads_to_bases(reads_N, pos)
        bases_T = reads_to_bases(reads_T, pos)
        
        counter_N = Counter(bases_N)
        counter_T = Counter(bases_T)
        
        a_N = counter_N[ref]
        b_N = counter_N[alt]
        d_N = a_N + b_N
        a_T = counter_T[ref]
        b_T = counter_T[alt]
        d_T = a_T + b_T
        
        if b_T <= 2 or b_N >= 2:
            continue
        
        minor_cn, major_cn = get_CN(pos, chrom_vcf, starts, ends, chroms, p_copy, m_copy)
        
        ID = chrom_vcf + ':' + str(pos)
        
        outcounts.write('\t'.join(map(str, [ID, a_T, b_T, 2, minor_cn, major_cn])) + '\n')
        
    
    
    normal_bam.close()
    tumor_bam.close()
    ref_genome_fasta.close()
    invcf.close()
    outcounts.close()
    

if __name__ == '__main__':
    main()
    
    
    
