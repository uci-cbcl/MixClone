#!/usr/bin/env python
import sys
import random
import numpy as np

centromere_starts = {}
centromere_ends = {}
heterochromatin_starts = {}
heterochromatin_ends = {}
q21_1_starts = {}
q21_1_ends = {}
excludes_starts = {}
excludes_ends = {}
centromere_starts['chr1'] = 121237000
centromere_ends['chr1'] = 123476957
heterochromatin_starts['chr1'] = 123476957
heterochromatin_ends['chr1'] = 141477000
q21_1_starts['chr1'] = 142400000
q21_1_ends['chr1'] = 148000000
excludes_starts['chr1'] = 121237000
excludes_ends['chr1'] = 148000000

chrom_start = 0

chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

chrom_lens = [247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
              158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
              114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
              63811651, 62435964, 46944323, 49691432]

subclone_prev_range = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

def get_genotype(max_copynumber):
    genotype = []
    
    for cn in range(0, max_copynumber+1):
        for B_num in range(0, cn+1):
            A_num = cn - B_num
            if A_num == 0 and B_num == 0:
                g_T = 'NULL'
            else:
                g_T = 'A'*A_num + 'B'*B_num
            
            genotype.append(g_T)
    
    return genotype


def get_genotype_edit_distance(g_x, g_y):
    A_num_x = g_x.count('A')
    B_num_x = g_x.count('B')
    A_num_y = g_y.count('A')
    B_num_y = g_y.count('B')
    
    edit_distance = np.abs(A_num_x - A_num_y) + np.abs(B_num_x - B_num_y)
    
    return edit_distance

def get_omega(G):
    tau = 100
    alpha = 0.5
    
    ed_genotype = get_genotype_edit_distance(G, 'AB')
    
    omega = np.power(alpha, ed_genotype)*tau
            
    return omega
    
def get_probs(genotype):
    omega = []
    
    for g in genotype:
        omega.append(get_omega(g))
    
    omega = np.array(omega)
    
    probs = omega*1.0/omega.sum()
    
    return probs

def get_P_M_copy(genotype):
    P_copy = []
    M_copy = []
    
    for g in genotype:
        P_copy.append(g.count('A'))
        M_copy.append(g.count('B'))
        
    return (P_copy, M_copy)



def main():
    outfile = open(sys.argv[1], 'w')
    chrom = sys.argv[2]
    max_copynumber = int(sys.argv[3])
    subclone_num = int(sys.argv[4])
    p_seg_num = int(sys.argv[5])
    q_seg_num = int(sys.argv[6])
    
    chrom_idx = chrom_list.index(chrom)
    p_start = chrom_start
    p_end = excludes_starts[chrom]
    q_start = excludes_ends[chrom]
    q_end = chrom_lens[chrom_idx]
    
    genotype = get_genotype(max_copynumber)
    genotype.remove('AB')
    genotype_probs = get_probs(genotype)
    P_copy, M_copy = get_P_M_copy(genotype)
    
    subclone_prevs = random.sample(subclone_prev_range, subclone_num)

    p_alpha = np.random.random_integers(1, p_seg_num, p_seg_num)
    q_alpha = np.random.random_integers(1, q_seg_num, q_seg_num)
    
    p_percent = p_alpha*1.0/p_alpha.sum()
    q_percent = q_alpha*1.0/q_alpha.sum()
    
    p_seg_lens = np.floor((p_end - p_start)*p_percent)
    q_seg_lens = np.floor((q_end - q_start)*q_percent)
    
    outfile.write('\t'.join(['#chrom', 'start', 'end', 'paternal_copy', 'maternal_copy', 'subclone_prev']) + '\n')
    print '\t'.join(['#chrom', 'start', 'end', 'paternal_copy', 'maternal_copy', 'subclone_prev'])
    
    base = p_start
    for i in range(0, p_seg_num):
        start = int(base)
        end = int(base + p_seg_lens[i])
        base = base + p_seg_lens[i]
        
        if i%2 == 0:
            p_copy = 1
            m_copy = 1
            subclone_idx = random.randint(0, subclone_num-1)
            subclone_prev = subclone_prevs[subclone_idx]
        else:
            subclone_idx = random.randint(0, subclone_num-1)
            genotype_idx = np.where(np.random.multinomial(1, genotype_probs) == 1)[0][0]
            p_copy = P_copy[genotype_idx]
            m_copy = M_copy[genotype_idx]
            subclone_prev = subclone_prevs[subclone_idx]
            
        outfile.write('\t'.join(map(str, [chrom, start, end, p_copy, m_copy, subclone_prev])) + '\n')
        print '\t'.join(map(str, [chrom, start, end, p_copy, m_copy, subclone_prev]))
        
    base = q_start
    for i in range(0, q_seg_num):
        start = int(base)
        end = int(base + q_seg_lens[i])
        base = base + q_seg_lens[i]
        
        if i%2 == 0:
            p_copy = 1
            m_copy = 1
            subclone_idx = random.randint(0, subclone_num-1)
            subclone_prev = subclone_prevs[subclone_idx]
        else:
            subclone_idx = random.randint(0, subclone_num-1)
            genotype_idx = np.where(np.random.multinomial(1, genotype_probs) == 1)[0][0]
            p_copy = P_copy[genotype_idx]
            m_copy = M_copy[genotype_idx]
            subclone_prev = subclone_prevs[subclone_idx]
            
        outfile.write('\t'.join(map(str, [chrom, start, end, p_copy, m_copy, subclone_prev])) + '\n')
        print '\t'.join(map(str, [chrom, start, end, p_copy, m_copy, subclone_prev]))
    
    outfile.close()
    
if __name__ == '__main__':
    main()