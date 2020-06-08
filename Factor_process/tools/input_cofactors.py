import numpy as np
import pandas as pd 
import gzip
from itertools import product
import os

import itertools as it
import time


###########################################################
########################################################### INPUT

from tools.compare_utilities import (
    get_available_muts, count_compare, deploy_count, pops_from_sim, check_availability, clean_empty
)

from tools.fasta_utilities import (
	geno_muts_v2, get_mutations, get_by_path,
	kmer_comp_index, kmer_mut_index
	)


##################################################
################################################## COUNT pop_kmers

def vcf_muts_matrix_v1(refseq,summary,start= 0,end= 0,ksize= 3,bases='ACGT', collapsed= True):
    ''' 
    Return matrix of mutation contexts by SNP in genotype array
    Each mutation is mapped to list of possible mutations as a binary vector.
    - v1 determines if alternative allele = reference allele in fasta. 
        if so, allele is switched, position idx is flagged. 
    '''
    
    mutations= get_mutations(bases= bases,ksize= ksize)
    kmers, kmer_idx= kmer_comp_index(mutations)
    
    mut_lib= kmer_mut_index(mutations)
    
    if end == 0:
        end= max(summary.POS)
    
    k5= int(ksize/2)
    k3= ksize - k5
    pos_mut= []
    flag_reverse= []
    flag_remove= []
    print('refseq len: {}'.format(len(refseq)))
    
    for x in range(summary.shape[0]):
        pos= int(summary.POS[x]) - 1
        if pos >=  start and pos <= end:
            kmer= refseq[pos-k5: pos + k3]
            if 'N' in kmer:
                flag_remove.append(x)
                continue

            if pos < k5 or (len(refseq) - pos) < k3:
                flag_remove.append(x)
                continue

            mut= kmer + summary.ALT[x]

            if kmer[1] == summary.ALT[x]:
                flag_reverse.append(x)
                mut= kmer+summary.REF[x]
            
            if len(mut) != 4: 
                print(kmer)
                print(summary.REF[x],summary.ALT[x])
                print(x,pos)
                print(len(refseq),summary.shape[0])
                if collapsed:
                    mut_array=np.zeros(len(kmer_idx))
                    pos_mut.append(mut_array)
                    continue
                else:
                    mut_array=np.zeros(len(mutations))
                    pos_mut.append(mut_array)
                    continue
            if collapsed:
                mut_index= kmers[mut]
                mut_array=np.zeros(len(kmer_idx))
            else:
                mut_index= get_by_path(mut_lib, list(mut))
                mut_array=np.zeros(len(mutations))
            
            mut_array[mut_index]= 1
            pos_mut.append(mut_array)
    
    pos_mut= np.array(pos_mut).T
    
    return pos_mut, flag_reverse, flag_remove


########################
######################## COUNT
########################

def lineAssign(ct,mut_idx,nmuts = 192):
    '''
    new_array= genotype vector. 
    mut_idx= indice of mutation type by SNP position. same order as new_array.
    nmuts= number of mutation types to be considered. 
    '''
    final= []
    for count_array in ct:
	    new_array= [0]*nmuts
	    for idx in range(len(count_array)):
	        nmut= mut_idx[idx]
	        new_array[nmut] += count_array[idx]

	    final.append(new_array)
    
    return np.array(final)


#######
def parse_PA(Window, pop_dict,frequency_range= [0,1],pop_tag= '_ss'):
    '''
    return dictionary with private alleles.
    '''
    PA_dict= {}

    pop_list= list(pop_dict.keys())
    print(Window.shape)
    pop_seg= {}
    pop_freqs= {}

    for pop in pop_list:
        t0= time.time()
        pop_ori= pop

        if pop_tag in pop:
            pop_ori= pop[len(pop_tag):].split('.')[0]

        klist= pop_dict[pop]
        pop_seg_ori= Window[klist,:]

        t1= time.time()
        pop_seg_ori= np.sum(pop_seg_ori,axis= 0)

        freqs= pop_seg_ori / len(klist)
        ## discount alleles outside freq range.
        in_out= (freqs <= frequency_range[0]) | (freqs >= frequency_range[1])

        pop_seg_ori= pop_seg_ori.reshape(1,Window.shape[1])
        #pop_seg_ori[:,in_out]= 0

        pop_seg_ori= pop_seg_ori > 0
        pop_seg_ori= np.array(pop_seg_ori,dtype= int).reshape(1,Window.shape[1])
        pop_seg[pop]= pop_seg_ori
        pop_freqs[pop]= in_out

    pop_array= [pop_seg[x] for x in pop_list]
    pop_array= np.array(pop_array)
    pop_sum= np.sum(pop_array,axis= 0)[0]

    PA_dict= {z: np.sum(pop_seg[z],axis= 0) for z in pop_list}
    PA_dict= {z: [x for x in range(len(g)) if g[x] == pop_sum[x] and pop_freqs[z][x] == False] for z,g in PA_dict.items()}

    return PA_dict

####### v1
def count_popKmers(Window, mut_matrix, mut_idx, pop_dict, single= True, frequency_range= [0,1],
									segregating= True, scale= 1, prop_gen_used= 1, return_private= False, return_seg= False,
									pop_tag= '_ss', row=32,col=3):
    '''
    Extract population mutation counts from _ind x kmer_ mutation matrix. 
    '''
    pop_counts= {}
    num_variants= {}
    pop_seg= {}
    PA_dict= {}

    pop_list= list(pop_dict.keys())

    for pop in pop_list:
        t0= time.time()
        pop_ori= pop

        if pop_tag in pop:
            pop_ori= pop[len(pop_tag):].split('.')[0]

        klist= sorted(pop_dict[pop])
        pop_gen= Window[klist,:]

        t1= time.time()

        if single: 
            pop_gen= np.sum(pop_gen,axis= 0)
            if segregating:
            	pop_gen= pop_gen > 0
                
            pop_gen= np.array(pop_gen,dtype= int).reshape(1,len(pop_gen))
        
        t2= time.time()
        if pop_gen.shape[0] == 1:
            pop_collapsed_mat= lineAssign(pop_gen,mut_idx,nmuts = mut_matrix.shape[0])
        else:
            pop_collapsed_mat= geno_muts_v2(pop_gen, mut_matrix)
        
        t3= time.time()

        tfetch= t1-t0
        tfilter= t2-t1
        tcount= t3-t2
        ttot= t3 - t0
        rate= ttot / (len(klist)/1000)
        print('#')
        print('w shape: {}'.format(Window.shape))
        print(pop)
        print('tfetch: {} s'.format(tfetch))
        print('t filter: {} s'.format(tfilter))
        print('N {} rate /1K : {}'.format(len(klist),rate))
        print('total {} s'.format(ttot))
        print('count {} s'.format(tcount / ttot))

        pop_seg[pop]= pop_gen

        pop_summed= np.sum(pop_collapsed_mat,axis= 0)
        t2= time.time()
        ######
        ######        
        pop_counts[pop]= pop_summed.reshape(row,col) * scale * prop_gen_used
        num_variants[pop]= np.sum(pop_collapsed_mat) * scale * prop_gen_used

    pop_summary= {
        'counts': pop_counts,
        'Nvars': num_variants,
        'sizes': {z:len(g) for z,g in pop_dict.items()}
    }

    if return_seg:
    	pop_summary['seg']= pop_seg

    return pop_summary, PA_dict

############################
#################### v2

def countkmers_cofactor(pop_gen,mut_matrix, mut_idx, pop_ori, single= True, frequency_range= [0,1],
    scale= 1, prop_gen_used= 1, return_private= False,PA= False):
    '''
    module to count_popKmers. this level allows to dissect populations.
    '''
    t0= time.time()
    if PA:
        freqs= np.sum(pop_gen,axis= 0) / pop_gen.shape[0]
        ## discount alleles outside freq range.
        in_out= (freqs <= frequency_range[0]) | (freqs >= frequency_range[1])
        pop_gen[:,in_out]= 0
    
    pop_seg_ori= np.sum(pop_gen,axis= 0) > 0
    pop_seg_ori= np.array(pop_seg_ori,dtype= int).reshape(1,len(pop_seg_ori))

    if single: 
        pop_gen= pop_seg_ori

    t1= time.time()
    if pop_gen.shape[0] == 1:
        pop_collapsed_mat= lineAssign(pop_gen[0],mut_idx,nmuts = mut_matrix.shape[0])
    else:
        pop_collapsed_mat= geno_muts_v2(pop_gen, mut_matrix)

    t2= time.time()
    print('#')
    print(pop_gen.shape)
    print(t1-t0)
    print(t2-t1)

    return pop_collapsed_mat, pop_seg_ori



def get_bins(array_len= 10, npacks= 1):
    '''
    get bins from array.
    '''
    t_range= np.arange(0,array_len,npacks,dtype= int)

    bins= []
    for idx in range(len(t_range)):
        if idx == len(t_range)-1:
            if t_range[idx] == array_len-1:
                bins[-1][1] = array_len-1
            else:
                bins.append([t_range[idx],array_len-1])
        else:
            bins.append([t_range[idx],t_range[idx+1]-1])

    return bins


def list_split(klist, npacks= 3):
    '''
    split list into smaller sections of at most len= npacks.
    '''
    bins= get_bins(array_len= len(klist), npacks= npacks)
    n_list= [klist[x[0]:(x[1]+1)] for x in bins]

    return n_list



def popCount_dissect(klist, Window, mut_matrix, pop_ori, single= True, frequency_range= [0,1],
                        scale= 1, prop_gen_used= 1, return_private= False,PA= {}, minN= 3e3):
    '''
    dissect population list to a smaller size to make operations faster.
    '''

    if len(klist) < minN:
        minN= len(klist)

    klist_diss= list_split(klist, npacks= minN)

    pop_gen= Window[klist_diss[0],:]

    pop_collapsed_mat, pop_seg_ori= countkmers_cofactor(pop_gen, mut_matrix, pop_ori, single= single, frequency_range= frequency_range,
                                            scale= scale, prop_gen_used= prop_gen_used, return_private= return_private,PA= PA)

    if len(klist_diss) > 1:

        for idx in range(1,len(klist_diss)):
            kset= klist_diss[idx]

            if not len(kset):
                continue

            pop_gen= Window[kset,:]

            pop_col, pop_seg= countkmers_cofactor(pop_gen,mut_matrix, pop_ori, single= single, frequency_range= frequency_range,
                                                    scale= scale, prop_gen_used= prop_gen_used, return_private= return_private,PA= PA)


            pop_collapsed_mat= np.concatenate((pop_collapsed_mat,pop_col),axis=0)
            pop_seg_ori+= pop_seg

        pop_seg_ori[pop_seg_ori > 0]= 1


    return pop_collapsed_mat, pop_seg_ori



########################################################### data management

def process_log(muted, sims_dir= ''):
    '''
    verify directories indicated in log exist in given directory.
    '''
    available= get_available_muts(muted)

    ### cleaning data set 
    ### i.e. accounting for aborted runs.
    available, miss_data= check_availability(available, dir_check=sims_dir)
    available, empty= clean_empty(available,str_format= '',dir_check= sims_dir,requested= ['.vcf.gz'])

    return available


def process_dir(sims_dir= ''):
    '''
    verify directories indicated in log exist in given directory.
    '''
    available= [name for name in os.listdir(sims_dir)]
    
    ### cleaning data set 
    ### i.e. accounting for aborted runs.
    available, miss_data= check_availability(available, dir_check=sims_dir)
    available, empty= clean_empty(available,str_format= '',dir_check= sims_dir,requested= ['.vcf.gz','fa.gz'])
    
    print('missing: {}, no vcf: {}'.format(len(miss_data),len(empty)))
    return available



def dict_write(new_dict,inds,outemp= 'ind_assignments{}.txt',dir_sim= '',tag= ''):
    '''
    cofactor to ind_assignment_scatter()
    '''
    inds= np.array(inds)
    new_array= [[(inds[x,0],v) for x in new_dict[v]] for v in new_dict.keys()]
    new_array= list(it.chain(*new_array))
    new_array= ['\t'.join(x) for x in new_array]
    
    out= dir_sim + outemp.format(tag)
    
    with open(out,'w') as f:
        f.write('\n'.join(new_array))


###################################
################################### Input without subsampling


def get_pop_dict(reference,dir_sim= '',indfile= 'ind_assignments.txt',haps_extract= False, return_inds= False):

    ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(ind_assignments,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    inds= [x for x in inds if x]
    pops= np.array(inds)[:,1]
    inds= np.array(inds)[:,0]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }

    total_N= len(inds)

    if haps_extract:
        pop_dict= {
            z: g + [x + total_N for x in g] for z,g in pop_dict.items()
        }
    
    if return_inds:
        return pop_dict, inds
    else:
        return pop_dict


####################################
#################################### Individual sub-sampling

def ind_assignment_scatter_v1(reference,dir_sim= '',indfile= 'ind_assignments.txt', haps_extract= False,
                          min_size= 80, samp= [5,20,10], stepup= "increment",outemp= 'ind_assignments{}.txt',
                          write_out= False, inds= [], pop_dict= {}, pop_sub= True):
    '''
    read ind assignments for a given window; 
    chose one population;
    subset that pop in some way.
    - v1: instead of writting new pop_assignment files, return them. 
    '''
    
    if not len(inds):
        ind_assignments= dir_sim + reference + '/' + indfile
        
        with open(ind_assignments,'r') as f:
            inds= f.readlines()
        
        inds= [x.split() for x in inds]
        pops= np.array(inds)[:,1]
        pop_dict= {
            z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
        }
        total_N= sum([len(x) for x in pop_dict.values()])
        
        if haps_extract:
            pop_dict= {
                z: g + [x + total_N for x in g] for z,g in pop_dict.items()
            }
        
    tag_list= []
    tag_dict= {}
    
    ## criterium of choice. chose only one pop.
    pop_avail= [x for x in pop_dict.keys() if len(pop_dict[x]) >= min_size]
    for pop_chose in pop_avail:
        
        N= len(pop_dict[pop_chose])
        pop_list= pop_dict[pop_chose]

        if pop_sub:
        	pop_list= list(range(len(pop_list)))

        if stepup== 'increment':
            timetable= np.linspace(2,samp[0],samp[1],dtype= int)
        else:
            timetable= np.arange(samp[0],N,samp[1],dtype= int) #np.linspace(samp[0],N,samp[1])

        print('timetable: len= {}; samp= {}'.format(len(timetable),samp))
        for each in timetable:  
            #each= int(each)
            for perm in range(samp[2]):
                tag= '_ss' + '.'.join([pop_chose,str(each),str(perm)])
                
                smaller= np.random.choice(pop_list,each,replace= False)
                new_pop= {
                    tag + '.s' + str(1):  smaller
                }
                
                #new_dict= {v:g for v,g in pop_dict.items() if v != pop_chose}
                #new_dict.update(new_pop)

                if write_out:
                    dict_write(new_pop,inds,outemp= outemp, dir_sim= dir_sim, tag= tag)
                else:
                    tag_dict[tag]= new_pop
                tag_list.append(tag)

    if write_out:
        return tag_list
    else: 
        return tag_list, tag_dict, pop_dict




def ind_assignment_dict(reference,pop_names, pop_lengths,dir_sim= '',indfile= 'ind_assignments.txt',
                          min_size= 80, samp= 20, stepup= "increment",outemp= 'ind_assignments{}.txt',write_out= False):
    '''
    read ind assignments for a given window; 
    chose one population;
    subset that pop in some way.
    - v1: instead of writting new pop_assignment files, return them. 
    '''
    
    ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(ind_assignments,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops)) if z in pop_names.keys()
    }

    pop_dict= {
        pop_names[z]: g for z,g in pop_dict.items()
    }
    
    tag_list= []
    tag_dict= {}
    
    ## criterium of choice. chose only one pop.
    pop_avail= list(pop_dict.keys())

    for perm in range(samp):
        
        new_dict= {}
        tag= '_ss' + str(perm)

        for pop_chose in pop_avail:

            N= pop_lengths[pop_chose]
            pop_list= pop_dict[pop_chose]

            smaller= np.random.choice(pop_list,N,replace= False)
            new_dict[pop_chose]= smaller


        tag_dict[tag]= new_dict
        tag_list.append(tag)

    if write_out:
        return tag_list
    else: 
        return tag_list, tag_dict, pop_dict




###############################
###############################


def read_diffs(tag,diff_dir= '',start= 0):
    '''
    read file of differences to ancestral sequence.
    '''

    filename= diff_dir + tag + '_diffs.txt.gz'

    with gzip.open(filename,'r') as f:
    	snps = f.readlines()
        
    snps= [x.decode() for x in snps[1:] if len(x)]
    snps= [x.split() for x in snps if 'SNP' in x]

    snps= {
    	str(int(x[0]) - start): [x[2],x[3]] for x in snps
    }

    return snps



def pop_dict_SFS(Window, pop_dict, ploidy= 2):
    '''
    read allele frequencies.
    '''

    
    pop_freqs= {}
    
    for pop in pop_dict.keys():
        pop_gen= Window[pop_dict[pop],:]
        pop_gen= np.sum(pop_gen,axis= 0)
        freq_dict= array_to_dict(pop_gen)
        #freq_dict= [(z,len(g)) for z,g in freq_dict.items()]

        pop_freqs[pop]= freq_dict
        
    return pop_freqs



from sklearn import decomposition
from sklearn.metrics import pairwise_distances


def pop_distances_PCA(Window, pop_dict, metric= 'euclidean', Npcs= 5, bandwidth= .2, dr= True):
    '''
    perform PCA on window numpy array, 
    parse transform coordinates by population indices in pop_dict. 
    calculate distances between centroids of populations by pairs.
    '''

    if dr: 
        pca = decomposition.PCA(n_components=Npcs)
        pca.fit(Window)
        
        data_pc= pca.transform(Window)
    else:
        data_pc= Window
    
    pop_data= {
        z: data_pc[g,:] for z,g in pop_dict.items()
    }
    pops= sorted(list(pop_dict.keys()))
    
    pop_combs= it.combinations(pops,2)
    pop_combs= list(pop_combs)
    
    comb_dict= {}
    
    for comb in pop_combs:
        pop1,pop2= comb
        centroids= {
            z: np.mean(g,axis= 0).reshape(1,-1) for z,g in pop_data.items()
        }

        distances= pairwise_distances(centroids[pop1],centroids[pop2],metric= metric)
        
        comb_dict[comb]= distances
    
    return comb_dict



