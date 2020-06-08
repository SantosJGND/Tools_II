import numpy as np
import itertools as it
from itertools import product
import os
import gzip

#from tools.plot_utilities import Population, frequency_breakdown, heatmap


def heatmap_mutation_labels():
    
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    ypos, ylabel = [], []

    mut_index = {}
    row, col = 0, 0

    labels= []

    for b2, d in [('A', 'T'), ('A', 'C'), ('A', 'G'),
                  ('C', 'T'), ('C', 'G'), ('C', 'A')]:

        for b1 in 'ACGT':
            row_lab= []
            col = 0
            ypos.append(row+0.5)
            if b1 == 'T' and b2 == 'C' and d == 'A':
                ylabel.append('5\'-'+b1)
            elif b1 == 'C':
                ylabel.append(b2+r'$\to$'+d+r'  '+b1)
            else:
                ylabel.append(b1)
            for b3 in 'ACGT':
                mut_index[(b1+b2+b3, d)] = (row, col)

                mut_index[(comp[b3]+comp[b2]+comp[b1], comp[d])] = (row, col)
                row_lab.append('_'.join([b1+b2+b3, d]))

                col += 1
            labels.append(row_lab)
            row += 1
    
    return labels



def get_available_muts(muted_log):
    ''' read log of mutation counts '''
    
    with open(muted_log,'r') as fp:
        available= fp.readlines()
    
    available= [x.strip() for x in available]
    
    return available


def pops_from_sim(sim,sims_dir= './mutation_counter/data/sims/',ind_file= "ind_assignments.txt",pop_set= True):
    '''read sim specific int to pop assignment, return pops.'''
    sim_dir= sims_dir + '{}/'.format(sim)
    ID_file= sim_dir + ind_file

    pops= []
    with open(ID_file,'r') as sample_id_lines:
        for line in sample_id_lines:
            line= str.encode(line)
            sample_id, population = line.split()[:2]
            pops.append(population.decode())
    
    if pop_set:
        return list(set(pops))
    else:    
        return pops


def count_compare(sim, frequency_range= [0,1], p_value= 1e-5, extract= 'pval', tag= '', muted_dir= './mutation_counter/data/mutation_count/',
                  sims_dir= './mutation_counter/data/sims/', exclude= False, ind_file= "ind_assignments.txt"):
    
    ''' perform pairwise population comparison of mutation counts for particular simulation'''
    pops= pops_from_sim(sim,sims_dir= sims_dir, ind_file= ind_file)

    ### change this 
    focus= pops[0]

    ## chromosome 
    chromosomes= [sim.split('.')[0].split('C')[1]]
    chromosome_groups = [chromosomes]

    ## get population pairs:
    population_pairs = [[(i), (i + 1) % len(pops)] for i in range(len(pops))] 
    population_pairs= [[pops[x] for x in y] for y in population_pairs]
    population_pairs= list(it.chain(*population_pairs))
    #print(population_pairs)
    pop_pair_names= zip(population_pairs[::2],
                               population_pairs[1::2])

    pop_pair_names= ['-'.join(list(x)) for x in pop_pair_names]

    population_pairs= zip(population_pairs[::2],
                               population_pairs[1::2])
    
    ### 
    chrom_pop= list(product(chromosome_groups,list(population_pairs)))

    heatmaps = [
        heatmap(
            chromosomes, population_pair, frequency_range, exclude, 
            p_value, sim, muted_dir, tag= tag, output= extract
        ) for chromosomes, population_pair in chrom_pop
    ]

    ratio_grids, significant_indices = zip(*heatmaps)
    
    return ratio_grids, significant_indices




def deploy_count(available, frequency_range= [0,1], p_value= 1e-5, extract= 'pval',muted_dir= './mutation_counter/data/mutation_count/',
                  sims_dir= './mutation_counter/data/sims/', tag= '', ind_file= "ind_assignments.txt"):
    
    ''' deploy count_compare() across simulations read from. '''
    data= {}
    
    for sim in available:
        
        ratio_grids, significant_indices= count_compare(sim, frequency_range= frequency_range, p_value= p_value, tag= tag,
                                               muted_dir= muted_dir, sims_dir= sims_dir, ind_file= ind_file)
        
        data[sim] ={
            'grids':ratio_grids,
            'sigs': significant_indices
        }
        
    return data



#####



def check_availability(available,str_format= '',dir_check= ''):
    '''check if names in list exist as directories somewhere, format possible.'''
    
    t= str_format.count('{}')
    
    if dir_check=='':
        dir_check= os.getcwd()
    
    if str_format:
        somlist= [str_format.format(*[x]*t) for x in available]
    else:
        somlist= list(available)
        
    dirs= [x[0] for x in os.walk(dir_check)]
    dirs= [x.split('/')[-1] for x in dirs]
    
    revised= [x for x in range(len(somlist)) if  somlist[x] in dirs]
    
    available= [available[x] for x in revised]
    missing= [available[x] for x in range(len(available)) if x not in revised]
    
    return available, missing


import subprocess

def gz_pattern(filename,patt= b'#CHROM'):
    '''
    check if pattern present in file.
    '''

    d= 0

    with gzip.open(filename,'rb') as fp:
        while d==0:
            nline= fp.readline()
            if patt in nline:
                d += 1

    return d


def clean_empty(available,str_format= '',dir_check= '',requested= ['.vcf.gz']):
    ''' check tag name specified directories for the presence of tag ID in files'''
    
    t= str_format.count('{}')
    
    if dir_check=='':
        dir_check= os.getcwd()
    
    if str_format:
        somlist= [str_format.format(*[x]*t) for x in available]
    else:
        somlist= list(available)
    
    av= []
    miss= []
    for idx in range(len(somlist)):
        sim= somlist[idx]
        ori= available[idx]
        dirhere=dir_check + sim + '/'
        files= [x[-1] for x in os.walk(dirhere)][0]
        
        present= []

        for fi in files:

            here= len([y for y in requested if y in fi])
            if here:
                if '.vcf.gz' in fi:
                    file = dirhere + fi
                    
                    ct= gz_pattern(file,patt= b"#CHROM")
                    if not ct:
                        continue

                present.append(fi)

        if len(present) >= len(requested):
            av.append(ori)
        else: miss.append(ori)
    
    return av, miss
        


def gzip_request(dir_check= '',str_format= '',requested= ['.vcf'], func= 'gzip'):
    ''' check tag name specified directories for the presence of tag ID in files'''
    
    t= str_format.count('{}')
    
    if dir_check=='':
        dir_check= os.getcwd()
    
    available= [name for name in os.listdir(dir_check)]
    
    if str_format:
        somlist= [str_format.format(*[x]*t) for x in available]
    else:
        somlist= list(available)
    
    av= []
    miss= []
    
    for idx in range(len(somlist)):
        sim= somlist[idx]
        ori= available[idx]
        dirhere=dir_check + sim + '/'
        files= [x[-1] for x in os.walk(dirhere)][0]
        
        present= [x for x in files if len([y for y in requested if y == x[-len(y):]]) > 0]
        
        if len(present):
            for fil in present:
                filename= dirhere + fil
                os.system('{} {}'.format(func,filename))
                miss.append(filename)
    
    print('{} files processed'.format(len(miss)))
    return miss

