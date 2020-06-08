
import re
import numpy as np
import pandas as pd 
import gzip
from itertools import product
import os

import itertools as it

import time

import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


from tools.fasta_utilities import (
	geno_muts_v2
	)


from tools.input_cofactors import (
		count_popKmers,
        vcf_muts_matrix_v1, process_log, process_dir, dict_write, 
		read_diffs, get_pop_dict, ind_assignment_scatter_v1, ind_assignment_dict,
		pop_distances_PCA, pop_dict_SFS
	)

#########################################
######################################### INPUT

def VCF_read_filter(sim, sim_dir= './',chrom= '1',haps_extract= False, scale_genSize= False,
    collapsed= True,min_size= 5, samp= [30,30,5], stepup= 'increment', outemp= './', ploidy= 2,
    indfile= 'ind_assignments.txt', fasta_template= 'chr{}_{}.fa.gz', diffs= False, bases= 'ACGT', ksize= 3):
    '''
    read vcf file. filter snps. 
    '''
    vcf_dir= sim_dir + sim + '/'
    vcf_file= vcf_dir + sim + '_' + 'chr' + chrom + '.vcf.gz'
    
    print(sim)

    #### read vcf
    genotype, summary, Names= read_vcf_allel(vcf_file,haps_extract= haps_extract)   
    
    print(genotype.shape)
    print(summary.shape)

    if len(genotype) == 0:
        return {}, {}, {}, {}, {}, {}
        
    ## read fasta
    fasta_file= vcf_dir + fasta_template.format(chrom,sim)

    with gzip.open(fasta_file,'r') as f:
        lines= f.readlines()
        lines= [x.decode() for x in lines]

    refseq= lines[1].strip()

    #### scale mutation counts. 
    L= len(refseq)
    scale= 1
    if scale_genSize:
        if genome_size > 1:
            scale= int(genome_size/L)

    ### subset window (for application to full data sets) if necessary.
    ### Not used here, wstart and wend set to min and max positions respectively.
    positions= [int(x) for x in summary.POS]
    wstart= int(min(positions))-1
    wend= int(max(positions))
    
    Wlen= wend - wstart
    
    genotype_parse= [x for x in range(summary.shape[0]) if int(summary.POS[x])-1 >= wstart and int(summary.POS[x])-1 <= wend]
    Window= genotype[:,genotype_parse]
    subset_summary= summary.loc[genotype_parse,:].reset_index()
    
    ### get mutation-type by SNP matrix, 
    ### filter SNPs if necessary (see vcf_muts_matrix_v1)
    t0= time.time()
    mut_matrix, flag_reverse, flag_remove= vcf_muts_matrix_v1(refseq,subset_summary,start= wstart,end= wend,ksize= ksize,
        bases=bases, collapsed= collapsed)
    
    print('mut_matrix shape: {}'.format(mut_matrix.shape))
    retain= [x for x in range(Window.shape[1]) if x not in flag_remove]
    Window= Window[:,retain]
    subset_summary= subset_summary.loc[retain,:].reset_index()

    t1= time.time()
    time_mut= t1 - t0

    if diffs:
    	sim_start= sim.split('.')[-1]
    	diff_snps= read_diffs(sim,diff_dir= vcf_dir, start= int(sim_start))

    	summary_diff= [x for x in range(subset_summary.shape[0]) if subset_summary.POS[x] in diff_snps.keys()]

    	flag_reverse.extend(summary_diff)
    	flag_reverse= list(set(flag_reverse))
    
    
    if flag_reverse:
        Window[:,flag_reverse]= ploidy - Window[:,flag_reverse]
    #

    return Window, mut_matrix, scale




import allel

def read_vcf_allel(file_vcf,haps_extract= False,calldata= 'calldata/GT'):
    '''
    Use scikit allel to read vcf file. Organise variant information into summary pandas df. 
    '''
    geno1= []

    vcf_ori= allel.read_vcf(file_vcf)

    if not vcf_ori:
        print('file:')
        print(file_vcf)
        print('is empty.')

        return {}, {}, {}

    ### get genotype array
    geno= vcf_ori[calldata]
    
    mult_alt= []
    indel= []
    single= []

    ## Filter SNPs. append to single list what to 
    for idx in range(geno.shape[0]):
        ## eliminate +1 segregating mutations.
        if vcf_ori['variants/ALT'][idx][1]:
            gen_t= geno[idx]
            gen_t[gen_t > 1] = 0
            geno[idx]= gen_t
            ## or just jump them
            indel.append(idx)

        elif len(vcf_ori['variants/REF'][idx]) != 1 or len(vcf_ori['variants/ALT'][idx][0]) != 1:
            indel.append(idx)
        else:
            single.append(idx)

    if haps_extract:
        geno1= geno[:,:,0].T
        geno= geno[:,:,1].T
        geno= np.concatenate((geno,geno1),axis= 0)
    else:
        geno= allel.GenotypeArray(geno)
        geno= geno.to_n_alt().T

    ## setup summary

    column_names= ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']

    alts= [vcf_ori['variants/ALT'][x][0] for x in range(geno.shape[1])]
    PASS= [['.','PASS'][int(vcf_ori['variants/FILTER_PASS'][x])] for x in range(geno.shape[1])]

    summary= [
        vcf_ori['variants/CHROM'],
        vcf_ori['variants/POS'],
        vcf_ori['variants/ID'],
        vcf_ori['variants/REF'],
        alts,
        vcf_ori['variants/QUAL'],
        PASS,

    ]

    summary= np.array(summary).T

    if len(indel):
        #
        geno= geno[:,single]
        if len(geno1):
            geno1= geno1[:,single]
        summary= summary[single,:]

    summary= pd.DataFrame(summary,columns= column_names)
    
    return geno, summary, vcf_ori['samples']


##################################################
##################################################




def MC_sample_matrix_v1(min_size= 80, samp= [5,20,10], stepup= "increment", diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',chrom_idx= 0, prop_gen_used= 1,
                    count_dir= './count/', dir_launch= '..',main_dir= './', sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/', segregating= False, scale_genSize= False,
                    outlog= 'indy.log', row= 24,col= 4, single= True, exclude= False, print_summ= False, sample_sim= 0,collapsed= True,bases= 'ACGT',ksize= 3,ploidy= 2, freq_extract= False, sim_del= 'C',
                    genome_size= 1,haps_extract= False, return_private= True):
    '''
    launch mutation counter pipeline on population assignments.
    Use matrix multiplication to extract counts. 
    - v1 relies on count_popKmers() function to count mutations per pop. allows freq. filter and single mutaiton count.  
    '''

    ti= time.time()
    sims= process_dir(sims_dir= sim_dir)
    print('available {}'.format(len(sims)))

    tags= []
    sim_extend= []
    chroms= []
    
    data_kmer= {}
    data_freqs= {}
    
    if sample_sim == 0:
        sample_sim= len(sims)
    
    print('sample {}'.format(sample_sim))
    sim_sub= np.random.choice(sims,sample_sim,replace= False)
    
    for sim in sim_sub:
        
        ## chromosome
        chrom= sim.split('.')[chrom_idx].split(sim_del)[-1].strip('chr')

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf
        t0= time.time()
        Window, mut_matrix, scale= VCF_read_filter(sim, sim_dir= sim_dir,chrom= chrom,haps_extract= haps_extract, scale_genSize= scale_genSize,
            collapsed= collapsed,min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,
            indfile= indfile,diffs= diffs,bases= bases, ksize= ksize, ploidy= ploidy)

        tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,dir_sim= sim_dir, haps_extract= haps_extract,
                      min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,indfile= indfile)
        
        total_inds= sum([len(x) for x in pop_dict.values()])
        t1= time.time()
        read_time= t1- t0
        if not len(Window) or Window.shape[0] < total_inds:
            continue

        ## counts for no tag sim:
        s0= time.time()
        pop_summary, PA_dict= count_popKmers(Window, mut_matrix, pop_dict, single= single, prop_gen_used= prop_gen_used,
                                  frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                  return_private= return_private)

        data_kmer[sim]= pop_summary

        if return_private: 
            pop_summary, dummy= count_popKmers(Window, mut_matrix, pop_dict, single= single, prop_gen_used= prop_gen_used,
                                      frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                      PA= PA_dict)
            data_kmer[sim]= pop_summary


        if freq_extract:
            pop_freqs= pop_dict_SFS(Window,pop_dict)
            data_freqs[sim]= pop_freqs
        
        t1= time.time()
        count_time= t1- t0
        
        if len(tag_list):
            ###
            sim_extend.append(sim)
            tags.append('')
            chroms.append(chrom)
            ###
            
            for idx in range(len(tag_list)):
                
                sim_extend.extend([sim]*len(tag_list))
                tags.extend(tag_list)
                chroms.extend([chrom]*len(tag_list))
                
                ##
                tag= tag_list[idx]
                ind_file= outemp.format(tags[idx])
                new_sim= sim + tag

                pop_dict= tag_dict[tag]
                
                pop_summary, dummy= count_popKmers(pop_summary['array'], mut_matrix, pop_dict, single= single, prop_gen_used= prop_gen_used,
                                  frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                  PA= PA_dict,counted= True)

                data_kmer[new_sim]= pop_summary

                if freq_extract:
                    pop_freqs= pop_dict_SFS(Window,pop_dict)
                    data_freqs[new_sim]= pop_freqs
                

        if print_summ:
            print('mut_matrix time: {} s'.format(time_mut / 60))
            print('count time: {} s'.format(count_time / 60))
            print('est total count time: {} s'.format(count_time*len(tag_list) / 60))
            print('replicates: {}'.format(len(tag_list)))
            print('read time: {} s'.format(read_time / 60))

    tf= time.time()
    time_elapsed= tf - ti
    
    print('time elapsed: {}s'.format(time_elapsed))

    return data_kmer, data_freqs



#####################
#####################

def MC_sample_matrix_simple(min_size= 80, samp= [5,20,10], stepup= "increment", diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',
                    count_dir= './count/', dir_launch= '..',main_dir= './', sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/', segregating= False,
                    outlog= 'indy.log', row= 24,col= 4, single= False, exclude= False, print_summ= False, sample_sim= 0,collapsed= True,bases= 'ACGT',ksize= 3,ploidy= 2, freq_extract= False,
                    distances= 'PCA',prop_gen_used= 1, scale_genSize= False, return_private= True,haps_extract= True):
    '''
    launch mutation counter pipeline on manipulated population assignments.
    Use matrix multiplication to extract counts. 
    - v1 relies on count_popKmers() function to count mutations per pop. allows freq. filter and single mutaiton count.  
    '''
    
    ti= time.time()
    sims= process_dir(sims_dir= sim_dir)
    print('available {}'.format(len(sims)))

    tags= []
    sim_extend= []
    chroms= []
    
    data_kmer= {}
    data_freqs= {}
    #
    if sample_sim == 0:
        sample_sim= len(sims)

    print('sample {}'.format(sample_sim))
    sim_sub= np.random.choice(sims,sample_sim,replace= False)
    
    for sim in sim_sub:

        ## chromosome
        chrom= sim.split('.')[0].split('C')[-1].strip('chr')
        chromosomes= [sim.split('.')[0].split('C')[1]]
        chromosome_groups = [chromosomes]

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf
        t0= time.time()
        Window, mut_matrix, scale= VCF_read_filter(sim, sim_dir= sim_dir,chrom= chrom,haps_extract= haps_extract, scale_genSize= scale_genSize,
            collapsed= collapsed,min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,
            indfile= indfile,diffs= diffs,bases= bases, ksize= ksize, ploidy= ploidy)

        tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,dir_sim= sim_dir, haps_extract= haps_extract,
                          min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,indfile= indfile)

        total_inds= sum([len(x) for x in pop_dict.values()])
        if not len(Window) or Window.shape[0] < total_inds:
            continue
        ## counts for no tag sim:
        s0= time.time()

        pop_summary, PA_dict= count_popKmers(Window, mut_matrix, pop_dict, single= single, prop_gen_used= prop_gen_used,
                                  frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                  return_private= return_private)

        data_kmer[sim]= pop_summary

        if return_private: 
            pop_summary, dummy= count_popKmers(Window, mut_matrix, pop_dict, single= single, prop_gen_used= prop_gen_used,
                                      frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                      PA= PA_dict)

            data_kmer[sim]= pop_summary

        if freq_extract:
            pop_freqs= pop_dict_SFS(Window,pop_dict)
            data_freqs[sim]= pop_freqs

        if distances:
            data_kmer[sim]['pairDist']= pop_distances_PCA(Window,pop_dict)
    
    return data_kmer, data_freqs



#######################################
#######################################



##################################################
################################################## 
#### READ
## read and sample from windows according to given dictionary. 


def MC_sample_matrix_dict(pop_names, pop_lengths,min_size= 80, samp= [5,20,10], stepup= "increment", diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',chrom_idx= 0,
                    count_dir= './count/', dir_launch= '..',main_dir= './', sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/', segregating= False, genome_size= 1,
                    outlog= 'indy.log', row= 24,col= 4, single= True, exclude= False, print_summ= False, sample_sim= 0,collapsed= True,bases= 'ACGT',ksize= 3,ploidy= 2, freq_extract= False, sim_del= 'C',
                         distances= 'PCA', Lsteps= 1,scale_genSize= False,prop_gen_used= 1,return_private= True):
    '''
    launch mutation counter pipeline on manipulated population assignments.
    Use matrix multiplication to extract counts. 
    - v1 relies on count_popKmers() function to count mutations per pop. allows freq. filter and single mutaiton count.  
    '''
    
    ti= time.time()
    sims= process_dir(sims_dir= sim_dir)
    print('available {}'.format(len(sims)))

    tags= []
    sim_extend= []
    chroms= []
    
    data_kmer= {}
    data_freqs= {}
    #sim_sample= np.random.choice(sims,8,replace= False)
    if sample_sim == 0:
        sample_sim= len(sims)
    
    print('sample {}'.format(sample_sim))
    sim_sub= np.random.choice(sims,sample_sim,replace= False)
    
    for sim in sim_sub:
        
        ## chromosome
        chrom= sim.split('.')[chrom_idx].split(sim_del)[-1].strip('chr')

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf
        t0= time.time()
        Window, mut_matrix, scale= VCF_read_filter(
        	sim, sim_dir= sim_dir,chrom= chrom,haps_extract= haps_extract, scale_genSize= scale_genSize,
            collapsed= collapsed,min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,
            indfile= indfile, diffs= diffs,bases= bases, ksize= ksize, ploidy= ploidy
            )

        tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,dir_sim= sim_dir, haps_extract= haps_extract,
                          min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,indfile= indfile)

        total_inds= sum([len(x) for x in pop_dict.values()])
        if not len(Window) or Window.shape[0] < total_inds:
            continue
        ## counts for no tag sim:
        s0= time.time()
        t1= time.time()
        count_time= t1- t0
        
        if len(tag_list):
            ###
            sim_extend.extend([sim]*len(tag_list))
            chroms.extend([chrom]*len(tag_list))
            ###
            Window_lengths= np.linspace(1,len(refseq),Lsteps,dtype= int)
            ###
            
            for idx in range(len(tag_list)):

                seq_idx= 0
                present_state= 0

                for snp_n in Window_lengths:
                    if snp_n < 10:
                        lrange= list(range(Window.shape[1]))
                        tag_l= 'full'
                    else:
                        while present_state < snp_n:
                            
                            if seq_idx >= (subset_summary.shape[0]-1):
                                present_state= len(refseq)
                                seq_idx= subset_summary.shape[0]-1
                            else:
                                present_state= subset_summary['POS'][seq_idx]
                                seq_idx += 1

                        lrange= list(range(seq_idx))
                        tag_l= str(snp_n * scale)
                    #
                    tag_here= tag_list[idx] + '.' + tag_l
                    tags.append(tag_here)
                    ##
                    tag= tag_list[idx]
                    #
                    new_sim= sim + tag_here

                    pop_dict= tag_dict[tag]
                    
                    #########
                    #########
                    pop_summary, PA_dict= count_popKmers(Window[:,lrange], mut_matrix[:,lrange], pop_dict, single= single, prop_gen_used= prop_gen_used,
                                              frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                              return_private= return_private)

                    data_kmer[new_sim]= pop_summary

                    if return_private: 
                        pop_summary, dummy= count_popKmers(Window[:,lrange], mut_matrix[:,lrange], pop_dict, single= single, prop_gen_used= prop_gen_used,
                                                  frequency_range= frequency_range,row=row,col=col,segregating= segregating,scale= scale,
                                                  PA= PA_dict)
                        data_kmer[new_sim]= pop_summary


                    if freq_extract:
                        pop_freqs= pop_dict_SFS(Window,pop_dict)
                        data_freqs[new_sim]= pop_freqs

                    if distances:
                        data_kmer[new_sim]['pairDist']= pop_distances_PCA(Window,pop_dict)

        if print_summ:
            print('mut_matrix time: {} s'.format(time_mut / 60))
            print('count time: {} s'.format(count_time / 60))
            print('est total count time: {} s'.format(count_time*len(tag_list) / 60))
            print('replicates: {}'.format(len(tag_list)))
            print('read time: {} s'.format(read_time / 60))

    tf= time.time()
    time_elapsed= tf - ti
    
    print('time elapsed: {}s'.format(time_elapsed))
    
    return data_kmer, data_freqs

