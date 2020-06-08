import re
import numpy as np
import pandas as pd 
import gzip
import itertools as it
from itertools import product
import os

from tools.compare_utilities import (
    get_available_muts, count_compare, deploy_count, pops_from_sim, check_availability, 
    clean_empty
)

from tools.fasta_utilities import (
    get_fasta_prop, get_mutations
    )

from tools.mcounter_cofactors import (
	set_SSD, get_pop_dict, get_chrom_sizes, 
	heatmap_v2, heatmap_v3, 
	run_stats, run_seg_stats
	)


import time
import scipy
import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


def md_SubSampVar(data, tag_ref= '_ss',row= 24,col= 4,bsep= 'C'):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii. for each population for each reference, organise first by popuation sizes (not proportions).
    iii. calculate pairwise differences between sets of counts are contiguous sample sizes. 
    '''
    avail= list(data.keys())
    print(avail.keys())
    
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail))]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }
    
    ref_pops= {}
    for idx in categ[0]:
        ref= avail[idx]
        for pop in data[ref]['counts'].keys():
            t= data[ref]['counts'][pop]
            t= t.reshape(1,np.prod(t.shape)) / np.sum(t)
            t= t[0]
            if pop in ref_pops.keys():
                ref_pops[pop].append(t)
            else:
                ref_pops[pop]=[t]
    
    for pop in ref_pops.keys():
        set1= ref_pops[pop]
        dists= set_SSD(set1,set1,same= True)
        ref_pops[pop]= dists
    
    #### population size diffs per population per simulation
    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}
    
    pops_vector= {}
    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        dat_size= [data[avail[av]]['sizes'][x] for x in dat]
        
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        dat_size= [dat_size[x] for x in range(len(dat))]
        dat_size= [round(x,3) for x in dat_size]
        
        for p in range(len(dat)):
            pops_vector[ref_pop[p]]= 0
            pop_asso[ref_sim][ref_pop[p]][dat_size[p]][avail[av]]= dat[p]
    
    pops_vector= list(pops_vector.keys())
    d= 0
    ### combine simulation combination and population size ranges.
    counts_dict= {}
    stats_dict= recursively_default_dict()
    for pop in pops_vector:
        counts_dict[pop]= {}
        for ref_sim in pop_asso.keys():
            
            batch= bsep.join(ref_sim.split(bsep)[:-1])
            available_sizes= sorted(list(pop_asso[ref_sim][pop].keys()))
            #
            size_counts= {}
            
            for si in available_sizes:
                t= [(v,g) for v,g in pop_asso[ref_sim][pop][si].items()]
                t= [data[v[0]]['counts'][v[1]] for v in t]
                t= [x.reshape(1,np.prod(x.shape)) / np.sum(x) for x in t if np.sum(x)]
                
                t= [x[0] for x in t]
                
                t= np.array(t)
                
                #
                if si in counts_dict[pop].keys():
                    counts_dict[pop][si].extend(t)
                else:
                    counts_dict[pop][si]=[t]
            
            counts_dict[pop]['sizes']= available_sizes
    
    for pop in counts_dict.keys():
        for size in counts_dict[pop].keys():
            set1= counts_dict[pop][size]
            
            dists_self= set_SSD(set1,set1,same= True)
            dists_ref= set_SSD(set1,ref_pops[pop],same= False)
            dists_ref= np.array(dists_ref).reshape(len(set1),len(ref_pops[pop]))
            print(dists_ref.shape)
            dists_ref= np.sum(np.array(dists_ref),axis= 1)

            stats_dict[pop][size]= {
                'self': dists_self,
                'ref': dists_ref
                }
    
    return stats_dict, counts_dict, ref_pops



############################################
############################################



def mcounter_deploy_v2(data,p_value= 1e-5, test_m= 'fisher', chi_total= False, ksize= 3, bases= 'ACGT',
                            frequency_range= [0,1], data_freqs= {}, extract= 'pval', row= 64,col= 3,
                            tag_ref= '_ss',collapsed=False,bsep= 'C',muted_dir= '', sims_dir= 'mutation_counter/data/sims/',
                            fasta_var= True):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii: for each pair of reference/subset populations, launch heatmapv2. return grid pvals or proportions,
    and proportion of mutations in subset population. allows for fisher or chi2 test for pval.
    - v2: compares sub pops to ref full pops other than its own; gets store of differences among refs.
    '''
    
    mutations= get_mutations(bases= bases,ksize= ksize)

    avail= list(data.keys())
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }
    
    fasta_ref_dict= {}
    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}
    ref_batch_dict= recursively_default_dict()

    reference_freqs= {}
    fasta_count_array= []
    
    for idx in categ[0]:
        ref= avail[idx]
        batch= bsep.join(ref.split(bsep)[:-1])
        ref_batch_dict[batch][ref]= ''
        ref_dir= sims_dir + ref + '/'
        fasta_kmer_prop= get_fasta_prop(ref,ref_dir,mutations,ksize= ksize,bases= bases,collapsed= collapsed)
        fasta_ref_dict[ref]= fasta_kmer_prop
        fasta_count_array.extend(fasta_kmer_prop)
        
    
    ### get fasta var
    fasta_count_array= np.array(fasta_count_array)
    print(fasta_count_array.shape)
    print('fasta array shape: {}'.format(fasta_count_array.shape))
    fasta_var= set_SSD(fasta_count_array,fasta_count_array,same= True)
    
    fasta_mtype_var= np.std(fasta_count_array,axis= 0)
    fasta_var_dict= {
        'dist': fasta_var,
        'type_sd': fasta_mtype_var
    }
    ###
    ref_batch_dict= {z:[] for z,g in ref_batch_dict.items()}
    
    if data_freqs:
        reference_freqs= {x: {} for x in pop_asso.keys()}
    
    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        for p in range(len(dat)):
            pop_asso[ref_sim][ref_pop[p]][avail[av]]= dat[p]

    d= 0
    count_data= recursively_default_dict()
    
    for ref in pop_asso.keys():
        batch= bsep.join(ref.split(bsep)[:-1])
        ref_count= fasta_ref_dict[ref]
        
        for pop in pop_asso[ref].keys():
            ref_pop_counts= data[ref]['counts'][pop]
            ref_pop_counts= ref_pop_counts / np.sum(ref_pop_counts)
            pop_shape= ref_pop_counts.shape
            ref_pop_counts= ref_pop_counts.reshape(1,np.prod(pop_shape))
            ref_pop_counts= (ref_pop_counts - ref_count / 3) / (1-ref_count)
            ref_pop_counts= ref_pop_counts.reshape(*pop_shape)
            ref_batch_dict[batch].append((pop,ref,ref_pop_counts))

            if data_freqs:
                reference_freqs[ref][pop]= data_freqs[ref][pop]
            
            for sub in pop_asso[ref][pop].keys():
                
                ref_pair= [(ref, pop),(sub, pop_asso[ref][pop][sub])]
                
                fasta_dict_local= {
                    ref: ref_count,
                    sub: fasta_ref_dict[ref]
                }
                
                count_data[d]= run_stats(ref,ref_pair,data,data_freqs= {},fasta_dict= fasta_dict_local,bsep=bsep,
                                                row= row,col= col,test_m= test_m, chi_total= chi_total)
                
                count_data[d]['pop']= pop
                count_data[d]['other']= []
                
                #for ref2 in pop_asso.keys():
                for pop2 in pop_asso[ref].keys():
                    if pop == pop2:
                        continue
                    if bsep.join(ref.split(bsep)[:-1]) != batch: 
                        continue
                    ##
                    pop_dict= {
                        ref: pop2,
                        sub: pop_asso[ref][pop][sub]
                    }
                    
                    fasta_dict_local= {
                        ref: ref_count,
                        sub: ref_count
                    } 
                    
                    ref_pair= [(ref, pop2),(sub, pop_asso[ref][pop][sub])]
                    
                    pair_stats= run_stats(ref,ref_pair,data,data_freqs= {},fasta_dict= fasta_dict_local,bsep=bsep,
                                                test_m= test_m,row= row,col= col,chi_total= chi_total)
                                            
                    count_data[d]['other'].append(pair_stats['diffs'])
                
                d += 1
    
    
    ####
    #### Reference counts, comparison by batch. 
    if len(fasta_var):
        return pop_asso, count_data, ref_batch_dict, reference_freqs, fasta_var_dict
    else:
        return pop_asso, count_data, ref_batch_dict, reference_freqs


##########################################################
########################################################## Binomial

def seg_comp_v2(data,p_value= 1e-5, test_m= 'fisher', Nbins= 20,stepup= 'increment', sampling= [],
                            frequency_range= [0,1], data_freqs= {}, extract= 'pval', chi_total= False,
                            muted_dir= '', tag_ref= '_ss',row= 64,col= 3):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii: for each pair of reference populations, launch heatmapv2. return grid pvals or proportions,
    and proportion of mutations in subset population. allows for fisher or chi2 test for pval.
    iii: v2 compares all possible combinations of sample sizes between population dyads.
    iv: v2 is no longer compatible with md_reference_comp. 
    '''
    
    Nmax= 1
    if stepup == 'increment':
        Nmax= sampling[0]
        #Nbins= sampling[0]
    
    bins= np.linspace(1,Nmax,Nbins)
    bins= np.round(bins,4)
    bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]
            
    avail= list(data.keys())
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }
        
    ### possible combinations per simulation.
    ref_combos= {}
       
    for idx in categ[0]:
        ref= avail[idx]
        ref_combs= list(data[ref]['seg'].keys())
        ref_combs= it.combinations(ref_combs,2)
        ref_combs= list(ref_combs)
        
        comb_dict= {
            x: {} for x in ref_combs
        }
        
        comb_stats= {}
        
        for pair in ref_combs:
            pop1, pop2= pair
            
            ref_pair= [(ref,pop1),(ref,pop2)]
            
            comb_stats[pair]= run_seg_stats(ref,ref_pair,data,data_freqs= data_freqs,row= row,col= col,
                                            chi_total= chi_total,return_counts= False)
        
        ref_combos[ref]= {
            'combs': comb_dict,
            'sizes': data[ref]['sizes'],
            'stats': comb_stats
        }
    
    #### population size diffs per population per simulation
    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}
    
    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        dat_size= [data[avail[av]]['sizes'][x] for x in dat]
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        dat_size= [dat_size[x] for x in range(len(dat))]

        if stepup != 'increment':
            dat_size= [dat_size[x] / data[avail[av]]['sizes'][ref_pop[x]] for x in range(len(dat))]
            dat_size= [round(x,3) for x in dat_size]
        for p in range(len(dat)):
            pop_asso[ref_sim][ref_pop[p]][dat_size[p]][avail[av]]= dat[p]
    
    d= 0
    ### combine simulation combination and population size ranges.
    
    for ref_sim in pop_asso.keys():
        batch= ref.split('C')[0]
        
        for combo in ref_combos[ref_sim]['stats'].keys():
            
            # sort availbale sizes for each population pair.
            # arrange them by bin.
            # deploy count comparisons for pairwise combinations of the two pops within each bin.
            
            pop1, pop2= combo
            
            available_sizes= {
                z: sorted(list(pop_asso[ref_sim][z].keys())) for z in combo
            }
            
            available_bins= {
                z: {
                    b: [x for x in available_sizes[z] if x > b[0] and x <= b[1]] for b in bins
                } for z in combo
            }
            
            ##
                        
            ij_queue= {}
            
            d= 0
            
            for i in available_bins[pop1]:
                for j in available_bins[pop2]:
                    
                    bend= (i,j)
                    #
                    ij_queue[bend]= []
                    ref_combos[ref_sim]['combs'][combo][bend]= {
                        'PA': {
                            v: [] for v in combo
                        },
                        'shared': {
                            v: [] for v in combo
                        }
                    }
                    
                    size_combs= [(x,y) for x in available_bins[pop1][i] for y in available_bins[pop2][j]]
                    #
                    for sc in size_combs:
                        s1, s2= sc
                        #
                        for sub1 in pop_asso[ref_sim][pop1][s1].keys():
                            for sub2 in pop_asso[ref_sim][pop2][s2].keys():
                                
                                ref_pair= [(sub1, pop_asso[ref_sim][pop1][s1][sub1]),(sub2, pop_asso[ref_sim][pop2][s2][sub2])]
                                
                                pop_count_dict= {
                                    z: [] for z in combo
                                }
                                
                                pop_counts= {
                                    g: data[g[0]]['seg'][g[1]] for g in ref_pair
                                }
                                
                                for idx in range(len(combo)):
                                    g= ref_pair[idx]
                                    
                                    pop_count_dict[combo[idx]].append(pop_counts[g][0])
                                    #ij_queue[bend].append(d)
                                    #d += 1
                    
                                ##
                                pop_count_dict= {
                                    z:np.array(g) for z,g in pop_count_dict.items()
                                }

                                PA_dict= {
                                    z: g + ref_combos[ref_sim]['stats'][combo]['PA'][(ref_sim,z)].T[:len(g.T)].T for z,g in pop_count_dict.items()
                                }

                                PA_dict= {
                                    z: np.array(g == 2,dtype= int) for z,g in PA_dict.items()
                                }

                                PA_dict= {z:np.sum(g) for z,g in PA_dict.items()}

                                shared_dict= {
                                    z: g + ref_combos[ref_sim]['stats'][combo]['shared'].T[:len(g.T)].T for z,g in pop_count_dict.items()
                                }

                                shared_dict= {
                                    z: np.array(g == 2,dtype= int) for z,g in shared_dict.items()
                                }

                                shared_dict= {z:np.sum(g) for z,g in shared_dict.items()}

                                for pop in combo:
                                    ref_combos[ref_sim]['combs'][combo][bend]['shared'][pop].append(shared_dict[pop])
                                    ref_combos[ref_sim]['combs'][combo][bend]['PA'][pop].append(PA_dict[pop])
            
    return pop_asso, ref_combos



