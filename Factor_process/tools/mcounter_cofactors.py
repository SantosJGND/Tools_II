

import numpy as np




def set_SSD(set1,set2,same= False):
    '''
    return sum of squared differences between every pair of vectors across two sets.
    '''
    dists= []
    
    for ind in set1:
        
        dist_vec= [(x - ind) for x in set2] #/ np.sum(indian + x)
        #dist_vec= [np.std(x) for x in dist_vec]
        dist_vec= [z**2 for z in dist_vec]
        dist_vec= [1 - np.sum(x) for x in dist_vec]
        dists.extend(dist_vec)
    
    if same:
        dists= np.array(dists)
        dists= dists.reshape(len(set1),len(set2))
        indx= np.triu_indices(len(set1),k=1)
        dists= dists[indx]
    
    return dists


def get_pop_dict(indfile= 'ind_assignments.txt',haps_extract= False):
    '''
    read and return ind pop assignment as dict.
    '''
    
    with open(indfile,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    inds= [x for x in inds if x]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }

    total_N= sum([len(x) for x in pop_dict.values()])

    if haps_extract:
        pop_dict= {
            z: g + [x + total_N for x in g] for z,g in pop_dict.items()
        }
    
    return pop_dict



def get_chrom_sizes(assembly= 'hg38', chrom_size_dir= '/home/jgarc235/Rhesus/chrom_sizes/', avoid= ['_','X','Y']):
    '''
    return dictionary of chrom sizes for given assembly
    '''

    filename= chrom_size_dir + assembly + '.chrom.sizes'

    with open(filename,'r') as fp: 
        lines= fp.readlines()
    
    lines= [x.split() for x in lines]
    lines= [x for x in lines if sum([int(y in x[0]) for y in avoid]) == 0]

    chrom_dict= {
        x[0].strip('chr'): int(x[1]) for x in lines
    }

    return chrom_dict


###### Heatmap variants - chi2
######
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact

def heatmap_v2(chromosomes,pop_counts, num_variants, population_dict,frequency_range= [0,1], exclude= False, 
                p_value= 5e-2, muted_dir= '',tag= '',output= 'pval',row= 64, col= 3, test= 'fisher',chi_total= False):
    '''
    pairwise comparison of count matrices. Chi2 applied cell-wise. 
    p-value or proportion - output argument. 
    - v2: count matrices are provided in pop_counts dictionary. 
    '''
    if test== 'fisher':
        from fisher import pvalue

    if exclude:
        files= read_exclude()
    else:
        files= {}
    
    refpop, pop = list(pop_counts.keys())

    ratio_grid = np.zeros((row, col))
    sig_x, sig_y = [], []
    
    for i in range(row):
        for j in range(col):
            if chi_total:
                gridsum= [np.sum(pop_counts[g]) for g in [pop,refpop]]
            else:
                gridsum= [np.sum(pop_counts[g][i]) for g in [pop,refpop]]

            chi_array= np.array([
                    [pop_counts[pop][i][j], gridsum[0] - pop_counts[pop][i][j]],
                    [pop_counts[refpop][i][j], gridsum[1] - pop_counts[refpop][i][j]]
                ],dtype= int)


            chi_0= np.sum(chi_array,axis= 1)
            chi_1= np.sum(chi_array,axis= 0)
            
            if chi_0[0] == 0 or chi_0[1] == 0:
                ratio_grid[i][j] = np.nan
                sig_x.append(j+0.5)
                sig_y.append(i+0.5)
            
            elif chi_1[0] == 0 or chi_1[1] == 0:
                ratio_grid[i][j] = 1
            
            else:
                ##
                if test == 'chi2':
                    _, this_pval, _, _ = chi2_contingency(
                        chi_array
                    )
                else:

                    p= pvalue(pop_counts[pop][i][j], np.sum(pop_counts[pop][i])  - pop_counts[pop][i][j],
                        pop_counts[refpop][i][j], np.sum(pop_counts[refpop][i]) - pop_counts[refpop][i][j])
                    this_pval= p.two_tail
                                    
                if output == 'pval':
                    ratio_grid[i][j] = this_pval
                else:
                    ratio_grid[i][j] = (pop_counts[pop][i][j] * num_variants[refpop] /
                                        (num_variants[pop] * pop_counts[refpop][i][j]))
                if this_pval < p_value:
                    sig_x.append(j+0.5)
                    sig_y.append(i+0.5)

    return ratio_grid, (sig_x, sig_y)



def heatmap_v3(chromosomes,pop_counts, num_variants, population_dict,frequency_range= [0,1], exclude= False, 
                p_value= 5e-2, muted_dir= '',tag= '',output= 'pval',row= 32, col= 3, test= 'fisher'):

    '''
    pairwise comparison of count matrices. Chi2 applied cell-wise. 
    p-value or proportion - output argument. 
    - v2: count matrices are provided in pop_counts dictionary. 
    - v3: performs chi-squared using mutation matrix rows instead of single cells.
    '''
    if exclude:
        files= read_exclude()
    else:
        files= {}
    
    refpop, pop = list(pop_counts.keys())

    ratio_grid = np.zeros((row, col))
    sig_x, sig_y = [], []
    
    for i in range(row):

        chi_array= np.array([
                pop_counts[pop][i],
                pop_counts[refpop][i]
            ],dtype= int)
        
        chi_0= np.sum(chi_array,axis= 1)
        chi_1= np.sum(chi_array,axis= 0)
        
        for idx in range(col):
            if chi_1[idx] == 0:
                chi_array[:,idx]= 1

        if chi_0[0] == 0 or chi_0[1] == 0:
            for j in range(col):
                ratio_grid[i][j] = np.nan
                sig_x.append(j+0.5)
                sig_y.append(i+0.5)

        elif chi_1[0] == 0 or chi_1[1] == 0:
            for j in range(col):
                ratio_grid[i][j] = 1

        else:
            ##
            if test == 'chi2':
                

                _, this_pval, _, _ = chi2_contingency(
                    chi_array
                )

            else:
                p= pvalue(chi_array[0][0],chi_array[0][1],
                    chi_array[1][0],chi_array[1][1])
                this_pval= p.two_tail

            for j in range(col):
                if output == 'pval':
                    ratio_grid[i][j] = this_pval
                else:
                    ratio_grid[i][j] = (pop_counts[pop][i][j] * num_variants[refpop] /
                                        (num_variants[pop] * pop_counts[refpop][i][j]))
                if this_pval < p_value:
                    sig_x.append(j+0.5)
                    sig_y.append(i+0.5)

    return ratio_grid, (sig_x, sig_y)


######
######
######



def run_stats(ref_sim,ref_pair,data,data_freqs= {},fasta_dict= {},bsep= 'C', test_m= 'fisher',
                    frequency_range= [0,1],exclude= False,p_value= 1e-2,row= 64,col= 3,muted_dir= '',chi_total= False):
    '''
    co-factor function to md counter comparisons, deploy heatmap and calculate kmer proportion differences 
    between pairs of population.
    - ref pair: list of tuples. can't be dictionary because of repeated pops / reference tags. 
    - correct mutation count frequency for mutation type freqs in fasta (fasta_dict). 

    '''
    batch= bsep.join(ref_sim.split(bsep)[:-1])
    sizes= [data[x[0]]['sizes'][x[1]] for x in ref_pair]
    #

    chromosomes= [ref_sim.split('.')[0].split('C')[1]]

    pop_counts= {
        g: data[g[0]]['counts'][g[1]] for g in ref_pair
    }

    num_variants= {
        g: data[g[0]]['Nvars'][g[1]] for g in ref_pair
    }

    num_variants= [num_variants[ref_pair[0]], num_variants[ref_pair[1]]]

    ratio_grid, sig_cells= heatmap_v2(chromosomes,pop_counts,num_variants,
                                      {},frequency_range= frequency_range, exclude= exclude, p_value= p_value,tag= '',
                                      test= test_m,output= 'pval',row= row,col= col,chi_total= chi_total)
    
    grid_props= pop_counts[ref_pair[1]] / pop_counts[ref_pair[0]]

    ## reshape arrays for pretty print
    shapes_dict= {z:s.shape for z,s in pop_counts.items()}
    pop_counts= {
        z: s.reshape(1,np.prod(s.shape)) for z,s in pop_counts.items()
    }
    ##

    ref_count= pop_counts[ref_pair[0]]

    pop_counts= {
        z: s / np.sum(s) for z,s in pop_counts.items()
    }
    #if fasta_dict:
    #    
    #    pop_counts= {
    #        z: (s - (fasta_dict[z[0]] / 3)) / (1-fasta_dict[z[0]]) for z,s in pop_counts.items()
    #    }
    
    ## reshape arrays for mutation type analysis later
    print(pop_counts[ref_pair[0]][:10])
    ##

    grid_diffs= pop_counts[ref_pair[0]] - pop_counts[ref_pair[1]]
    pop_ref_prop= pop_counts[ref_pair[0]]
    pop_ref_prop[pop_ref_prop == 0] = 1

    print(grid_diffs[0,:10])

    grid_diffs= grid_diffs / pop_ref_prop
    print(grid_diffs[0,:10])
    print('#')

    grid_diffs= grid_diffs.reshape(*shapes_dict[ref_pair[0]])
    pop_counts= {z: s.reshape(*shapes_dict[z]) for z,s in pop_counts.items()}

    comb_stats= {
        'grids': ratio_grid,
        'sigs': sig_cells,
        'sizes': sizes,
        'batch': batch,
        'diffs': grid_diffs,
        'props': grid_props, 
        'counts': pop_counts,
        'Nvar': num_variants
    }
    
    if data_freqs:
        comb_stats['freqs']= {
            x: data_freqs[x[0]][x[1]] for x in ref_pair
        }
    
    return comb_stats




def run_seg_stats(ref_sim,ref_pair,data,data_freqs= {},fasta_dict= {},bsep= 'C', test_m= 'fisher',
                    frequency_range= [0,1],exclude= False,p_value= 1e-2,row= 64,col= 3,muted_dir= '',
                  chi_total= False,return_counts= True):
    '''
    segregating variant stats. extracts shared and private alleles between ref_pair pops. 
    '''
    batch= ref_sim.split(bsep)[0]
    #

    pop_counts= {
        g: data[g[0]]['seg'][g[1]] for g in ref_pair
    }
    
    count_array= [pop_counts[g] for g in ref_pair]
    count_array= count_array[0] - count_array[1]
    count_array= count_array[0]
    
    num_priv= {
        ref_pair[z]: np.array(count_array == [1,-1][int(z)],dtype= int) for z in [0,1]
    }
    
    total= np.array(count_array == 0,dtype= int)
    
    
    if return_counts:
        total= np.sum(total)
        num_priv= {z: np.sum(g) for z,g in num_priv.items()}
    
    comb_stats= {
        'PA': num_priv,
        'shared': total
    }
        
    return comb_stats


