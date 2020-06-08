import numpy as np
from itertools import product
import os
from functools import reduce  # forward compatibility for Python 3
import operator
import gzip
import tempfile
import collections


def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)



def reference_sequence(chromosome_number,reference= 'rheMac10',dir_launch='..'):
    '''
    read single chromosome fasta.
    '''
    file_path_template = dir_launch+'chr{}_{}.fa.gz'.format(chromosome_number,reference)
    
    with gzip.open(file_path_template) as infile:
        lines = infile.readlines()
    
    infile.close()
    
    f = tempfile.TemporaryFile()
    
    for v in range(1,len(lines)):
        processed_line= lines[v].upper().strip(b'\n')

        f.write(processed_line)
    
    f.seek(os.SEEK_SET)
    result= f.read()
    
    f.close()
    
    return result

####
####
def get_mutations(bases= 'ACGT',ksize= 3):
    '''return list of possible kmer mutations'''
    
    mutations=[]
    
    base_set= [bases]*ksize

    for trimer in product(*base_set):
        for base in bases:
            if trimer[int(ksize / 2)] != base:
                mutations.append((''.join(trimer), base))
    
    return mutations

####
def get_by_path(root, items):
    """Access a nested object in root by item sequence."""
    return reduce(operator.getitem, items, root)

def set_by_path(root, items, value):
    """Set a value in a nested object in root by item sequence."""
    get_by_path(root, items[:-1])[items[-1]] = value

####
####

def kmer_dict_init(ksize= 3,bases='ATCG'):
    '''produce nested dictionary of nucs for a particular kmer size'''
    mut_lib= recursively_default_dict()

    base_set= [bases]*ksize

    for trimer in product(*base_set):
        set_by_path(mut_lib,list(trimer),0)
    
    return mut_lib


def fasta_get_freq(seq,start= 0,end= 0,step= 1,ksize=3,bases= 'ACGT'):
    '''return count of kmer across fasta region'''
    kmer_dict= kmer_dict_init(ksize= ksize,bases=bases)
    if end == 0:
        end= len(seq) - ksize
    
    for ki in range(start,end,step):
        kmer= seq[ki:ki+ksize]
        if 'N' in kmer:
            continue
        get_by_path(kmer_dict, kmer[:-1])[kmer[-1]] += 1
    
    return kmer_dict




def get_complement(kmer):
    '''Return complement of a given kmer'''
    complements= {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    
    comp= [complements[x] for x in kmer][::-1]
    return comp



def complement_dicts(bases= 'ACGT',ksize= 3):
    '''return dict of comp kmers + index dict to parse with'''
    comp_dict= {}
    comp_index= {}
    d= 0
    base_set= [bases] * ksize
    mers= product(*base_set)
    
    for kmer in mers:
        kmer= ''.join(kmer)
        if kmer not in comp_dict.keys():
            comp_dict[kmer]= d
            
            comp= get_complement(kmer)
            comp= ''.join(comp)
            comp_dict[comp]= d
            
            comp_index[d]= (kmer,comp)
            d += 1
    
    return comp_dict,comp_index



def fasta_get_freq(seq,start= 0,end= 0,step= 1,ksize=3,bases= 'ACGT'):
    '''return count of kmer across fasta region'''
    kmer_dict= kmer_dict_init(ksize= ksize,bases=bases)
    if end == 0:
        end= len(seq) - ksize

    for ki in range(start,end,step):
        kmer= seq[ki:ki+ksize]
        if 'N' in kmer:
            continue
        get_by_path(kmer_dict, kmer[:-1])[kmer[-1]] += 1

    return kmer_dict


def kmer_freq_balance(kmer_dict, mutations, fasta_len= 10000, bases= 'ACGT',ksize= 3,collapsed= False):
    '''return list of possible kmer mutations'''
    
    kmers, kmer_idx= kmer_comp_index(mutations)
    
    mutation_sum= [0] * int(len(mutations) / (1 + int(collapsed)))
    Nkmers= fasta_len - ksize

    for idx in range(len(mutations)):
        mut= mutations[idx]
        prop= get_by_path(kmer_dict,mut[0])
        prop= prop / Nkmers
        
        if collapsed:
            idx= kmers[''.join(mut)]
        
        mutation_sum[idx]= prop

    return np.array(mutation_sum).reshape(1,-1)




def get_fasta_prop(sim,sim_dir,mutations,ksize= 3,bases= 'ACGT',collapsed= False,anc_corr= '',bsep= 'C'):
    '''
    Count proportion of kmers along a fasta sequence. 
    Allow to collapse kmer counts by reverse complementarity. 
    '''
    
    chrom= sim.split('.')[0].split(bsep)[-1]
    
    fasta_file= sim_dir + 'chr{}_{}.fa.gz'.format(chrom,sim)

    with gzip.open(fasta_file,'r') as f:
        lines= f.readlines()
        lines= [x.decode() for x in lines]

    refseq= lines[1].strip()

    kmer_dict= fasta_get_freq(refseq,start= 0,end= 0,step= 1,ksize=ksize,bases= bases)
    
    ref_kmer_prop =kmer_freq_balance(kmer_dict,mutations,fasta_len= len(refseq),collapsed= collapsed)
    
    return ref_kmer_prop




###########
########### get individual kmer frequencies following mutation.

def collapse_freqs(kmer_dict,comp_index):
    '''return vector of collapsed counts by kmer'''
    counts= []
    
    for kdx in comp_index.keys():
        total= [get_by_path(kmer_dict, list(comp)) for comp in comp_index[kdx]]
        total= sum(total)
        counts.append(total)
        
    return counts
            

### vcf - fasta functions

def vcf_kmers(refseq, summary, comp_dict, ksize=3,start= 0,end= 0):
    '''
    get kmers in vcf.
    returns list of reference and alternative kmer indices in comp_dict.
    '''
    if end == 0:
        end= max(summary.POS)
        
    ###
    kmer_dict= {}
    kmer_list= []
    d= 0
    base_set= [bases] * ksize
    mers= product(*base_set)
    
    for kmer in mers:
        kmer= ''.join(kmer)
        kmer_dict[kmer]= d

    ###
    k5= int(ksize/2)
    k3= ksize - k5
    posKmer_ref= []
    posKmer_alt= []
    
    for x in range(summary.shape[0]):
        pos= int(summary.POS[x]) - 1
        if pos >=  (start + k5) and pos <= (end - k3):
            kmer= refseq[pos-k5: pos + k3]
            mut= kmer[:k5] + summary.ALT[x] + kmer[k3:]
            
            posKmer_ref.append(comp_dict[kmer])
            posKmer_alt.append(comp_dict[mut])
    
    return posKmer_ref, posKmer_alt



######################################
######################################


def geno_kmers(genotype, summary, refseq,comp_index,ksize= 3,bases= 'ATCG', start= 0, end= 0):
    '''get individual collapsed mutation arrays across data set.'''
    
    k5= int(ksize/2)
    k3= ksize - k5
    
    ind_dicts= {ind: fasta_get_freq(refseq,start= int(start),end= int(end),step= 1,ksize=ksize,bases= bases) for ind in range(genotype.shape[0])}
    collapsed= []
    
    for snp in range(genotype.shape[1]):
        alleles= genotype[:,snp]
        
        changes= [x for x in range(len(alleles)) if alleles[x] > 0]
        pos= int(summary.POS[snp]) - 1
        
        kmer= refseq[pos-k5: pos+k3]
        kmer_comp= get_complement(kmer)
        mut= kmer[:k5] + summary.ALT[snp] + kmer[k3:]
        mut_comp= get_complement(mut)
        
        for ind in changes:
            get_by_path(ind_dicts[ind], kmer[:-1])[kmer[-1]]-=  alleles[ind]
            get_by_path(ind_dicts[ind], mut[:-1])[mut[-1]]+= alleles[ind]
    
    for ind in range(genotype.shape[0]):
        collapsed_freqs= collapse_freqs(ind_dicts[ind],comp_index)
        collapsed.append(collapsed_freqs)
    
    collapsed= np.array(collapsed)
    collapsed= (collapsed.T/collapsed.sum(axis=1)).T
    
    return collapsed




############################################
############################################ Using dictionaries.


def kmer_mut_init(mutations, default= 0):
    '''produce nested dictionary of nucs for a particular kmer size'''
    
    mut_lib= recursively_default_dict()
    
    d= 0
    for mut in range(len(mutations)):
        trimer= mutations[mut]
        trimer= ''.join(trimer)
        get_by_path(mut_lib, trimer[:-1])[trimer[-1]]= 0
    
    return mut_lib


def kmer_mut_index(mutations):
    '''produce nested dictionary of nucs for a particular mutation list'''
    mut_lib= recursively_default_dict()
    
    for mut in range(len(mutations)):
        trimer= ''.join(mutations[mut])
        get_by_path(mut_lib, trimer[:-1])[trimer[-1]]= mut
    
    return mut_lib


def vcf_muts(refseq,summary,start= 0,end= 0,ksize= 3,bases='ATCG'):
    ''' return vector of mutation contexts by SNP in vcf. '''
    
    mut_lib= kmer_mut_index(mutations)
    
    if end == 0:
        end= max(summary.POS)
    
    k5= int(ksize/2)
    k3= ksize - k5
    pos_mut= []
    
    for x in range(summary.shape[0]):
        pos= int(summary.POS[x]) - 1
        if pos >=  start and pos <= end:
            kmer= refseq[pos-k5: pos + k3]
            mut= kmer + summary.ALT[x]
            mut_index= get_by_path(mut_lib, list(mut))
            
            pos_mut.append(mut_index)
            
    return pos_mut



def kmer_comp_index(mutations):
    ''' return nested dictionaries of kmer mutations w/ index'''
    kmers= {}
    kmer_idx= {}
    d= 0
    for kmer in mutations:

        comp= get_complement(kmer[0]) + get_complement(kmer[1])
        comp= ''.join(comp)
        kmer= ''.join(kmer)
        
        if comp in kmers.keys():
            idx= kmers[comp]
            kmers[kmer]= idx
            kmer_idx[idx].append(kmer)
        else:
            kmers[kmer]= len(kmer_idx)
            kmer_idx[len(kmer_idx)]= [kmer]

        d += 1
    
    return kmers, kmer_idx


def collapse_muts(kmer_dict,mutations, collapse= True):
    '''return vector of counts by kmer, optional collapse by complement'''
    counts= []
    
    if collapse:
        kmers, kmer_idx= kmer_comp_index(mutations)

        for kdx in kmer_idx.keys():
            count= [get_by_path(kmer_dict, list(comp)) for comp in kmer_idx[kdx]]
            count= sum(count)
            ##
            counts.append(count)
    
    else:
        for comp in mutations:
            comp= ''.join(comp)
            kdx= get_by_path(kmer_dict, list(comp))
            counts.append(kdx)
        
    return counts


def geno_muts_v1(geno_array, kmer_dict, vcf_muts_vector, mutations, 
                 bases= 'ACGT', ksize= 3, Wl= 0, collapse= True):
    ''' return mutation spectrum'''
    
    ind_dict_store= {
        i: kmer_mut_init(mutations) for i in range(geno_array.shape[0])
    }
    collapsed= []
    
    for snp in range(geno_array.shape[1]):
        
        alleles= geno_array[:,snp]
        changes= [x for x in range(len(alleles)) if alleles[x] > 0]
        
        minus= mutations[vcf_muts_vector[snp]]
        
        prop_kmer= get_by_path(kmer_dict, minus[0]) / (Wl - ksize + 1)
        minus= ''.join(minus)        
        
        for ind in changes:
            get_by_path(ind_dict_store[ind], minus[:-1])[minus[-1]]+= alleles[ind] #/ prop_kmer
    
    for ind in range(geno_array.shape[0]):
        collapsed_freqs= collapse_muts(ind_dict_store[ind],mutations,collapse= collapse)
        collapsed.append(collapsed_freqs)
    
    collapsed= np.array(collapsed)
    collapsed= (collapsed.T/collapsed.sum(axis=1)).T
    
    return collapsed

########################################################
######################################################## Using matrices


def vcf_muts_matrix(refseq,summary,start= 0,end= 0,ksize= 3,bases='ATCG', collapse= True):
    ''' 
    Return matrix of mutation contexts by SNP in genotype array
    Each mutation is mapped to list of possible mutations as a binary vector.
    '''
    
    mutations= get_mutations(bases= bases,ksize= ksize)
    kmers, kmer_idx= kmer_comp_index(mutations)
    
    mut_lib= kmer_mut_index(mutations)
    
    if end == 0:
        end= max(summary.POS)
    
    k5= int(ksize/2)
    k3= ksize - k5
    pos_mut= []
    
    for x in range(summary.shape[0]):
        pos= int(summary.POS[x]) - 1
        if pos >=  start and pos <= end:
            kmer= refseq[pos-k5: pos + k3]
            mut= kmer + summary.ALT[x]
            
            if collapse:
                mut_index= kmers[mut]
                mut_array=np.zeros(len(kmer_idx))
            else:
                mut_index= get_by_path(mut_lib, list(mut))
                mut_array=np.zeros(len(mutations))
            
            mut_array[mut_index]= 1
            pos_mut.append(mut_array)
    
    pos_mut= np.array(pos_mut).T
    
    return pos_mut



def geno_muts_v2(geno_array, vcf_muts_matrix, standardize= False):
    ''' 
    Return mutation spectrum using matrix multiplication.
    multiply mutation matrix by genotype array to obtain matrix of mutations by samples.
    '''
    
    collapsed= vcf_muts_matrix @ geno_array.T
    
    collapsed= np.array(collapsed).T

    if standardize:
        collapsed= (collapsed.T/collapsed.sum(axis=1)).T
    
    return collapsed
