
import numpy as np 

import os
import time

from tools.input_cofactors import (
	ind_assignment_scatter_v1, count_popKmers
	)


from tools.mcounter_cofactors import (
	set_SSD
	)


##################
def read_simCounts(simdb,tag_split= 'C',pop_tag= '_ss'):
	'''
	read file of individual mutation type counts. 
	first 3 columns= simID, pop, ind. 
	header= True.
	'''
	with open(simdb,'r') as fp:
		counts= fp.readlines()

	header= counts[0]
	muts= header[3:]
	counts= [x.strip().split('\t') for x in counts[1:]]
	counts= np.array(counts)
	pop_names= counts[:,1]
	for idx in range(len(pop_names)):
		pop= pop_names[idx]
		if pop_tag in pop:
			pop= pop[len(pop_tag):].split('.')[0]
			pop_names[idx]= pop 

	counts[:,1]= pop_names

	return counts, muts, header




def info_array_collect(db_dir, 
					row= 24,col= 4, tag_split= 'C', tag_sim= '_ss'):
	'''
	extract mutation count sub-samples from across simulation count dbs. 
	'''

	list_neigh= os.listdir(db_dir)
	lines= []
	print(list_neigh)

	for dbf in list_neigh[:6]:
		simdb= db_dir + dbf

		counts, muts, header= read_simCounts(simdb)

		batch= counts[0][0]
		print('base: {}'.format(batch))
		print(counts.shape)
		#print(counts[:10,:10])

		lines.append(counts)

	lines= np.concatenate(tuple(lines),axis= 0)

	info_array= lines[:,:3]

	counts= lines[:,3:]
	counts= np.array(counts,dtype= int)

	return info_array, muts, counts



##############
############## Managing.
############## Per simulation analysis. 
#### I. Cofactor functions.

def sim_countPrep(count_array, info_array, sim, sim_idx):
	'''
	return count proportion array and extract factor information from count matrix. 
	Make prettier possible w/ functional factor get.
	'''
	if len(sim_idx) == 1:
		return "","","","",""

	lines= count_array[sim_idx,:]
	#print(sim_idx)
	print(count_array.shape)
	ncounts= np.sum(lines,axis= 1)

	info_sim= info_array[sim_idx,:]

	pop_file= info_sim[:,1]
	size_file= np.array(info_sim[:,2],dtype= int)
	pop_list= list(set(pop_file))

	props= lines.T / ncounts
	props= props.T
	#print(sim)
	#print(props.shape)
	#print(pop_list)

	return props, pop_list, size_file, pop_file



def stats_condense(props, pop_list, size_file, pop_file):
	'''
	subset count data by population and size. 
	Assume that population with largest size is to be taken as reference. 
	get count proportions and differences to population and simulation specific reference data. 
	'''
	### pop and size index and size dict
	pop_idx_dict= {
	    pop: [x for x in range(len(pop_file)) if pop_file[x] == pop] for pop in pop_list
	}

	pop_size_dict= {
	    pop: [size_file[x] for x in pop_idx if pop_file[x] == pop] for pop,pop_idx in pop_idx_dict.items()
	}


	pop_size_idx_dict= {
	    pop: {
	        si: [x for x in range(len(pop_si)) if pop_si[x] == si] for si in list(set(pop_si)) 
	    } for pop,pop_si in pop_size_dict.items()
	}

	### ref dicts
	pop_refkeys= {
	    pop: max(list(g.keys())) for pop,g in pop_size_idx_dict.items()
	}

	### count and diffs arrays:

	pop_props= {
	    pop: props[g,:] for pop,g in pop_idx_dict.items()
	}

	ref_props= {
	    pop: pop_size_idx_dict[pop][pop_refkeys[pop]] for pop in pop_list
	}

	ref_props= {
	    pop: pop_props[pop][g,:] for pop,g in ref_props.items()
	}

	ref_props= {
	    pop: np.mean(g,axis= 0) for pop,g in ref_props.items()
	}

	pop_diffs= {
	    pop: (ref_props[pop] - g) / ref_props[pop] for pop,g in pop_props.items()
	}

	return pop_diffs, ref_props, pop_props, pop_refkeys, pop_size_idx_dict, pop_size_dict, pop_idx_dict


##########
########## Counts for specific purposes. 
##########

def sim_data(count_array,info_array,row= 48,col= 4):
	'''
	extract count differences, proportions and counts from mutation type count array.
	labels provided: sim, pop_labels and population size labels = columns 0,1 & 2 of info_array. 
	pop_list= list of populations to extract. must exist in pop_labels.
	'''

	sim_labels= list(info_array[:,0])
	sim_dict= {x:[] for x in list(set(sim_labels))}
	print(sim_dict.keys())
	for idx in range(len(sim_labels)):
		sim_dict[sim_labels[idx]]+= [idx]

	#####
	d= 0

	count_data= {}
	count_sims= list(sim_dict.keys())
	for sim,sim_idx in sim_dict.items():

		props, pop_list, size_file, pop_file= sim_countPrep(count_array, info_array, sim, sim_idx)

		if not len(props):
			continue

		pop_diffs, ref_props, pop_props, pop_refkeys, pop_size_idx_dict, pop_size_dict, pop_idx_dict= stats_condense(props, pop_list, size_file, pop_file)
		
		for pop,si_dict in pop_size_idx_dict.items():
		    for si,si_idx in si_dict.items():
		        for idx in si_idx:
		            local_array= count_array[pop_idx_dict[pop][idx],:]
		            if np.nanmin(pop_props[pop][idx]) < 0:
		            	continue
		            print('#')
		            print(pop)
		            print(si)
		            print(local_array[:20])
		            print(sum(local_array))
		            count_data[d]= {
		                'pop': pop,
		                'sizes': [pop_refkeys[pop],si],
		                'Nvar': [0,np.sum(local_array)],
		                'props': pop_props[pop][idx].reshape(row,col),
		                'diffs': pop_diffs[pop][idx].reshape(row,col)
		            }
		            d+=1

	return count_data


def sim_VarSub(count_array,info_array,row= 48,col= 4,si_max= 100):
	'''
	extract count differences, proportions and counts from mutation type count array.
	labels provided: sim, pop_labels and population size labels = columns 0,1 & 2 of info_array. 
	pop_list= list of populations to extract. must exist in pop_labels.
	'''

	sim_labels= list(info_array[:,0])
	sim_dict= {x:[] for x in list(set(sim_labels))}
	print(len(sim_dict.keys()))
	for idx in range(len(sim_labels)):
		sim_dict[sim_labels[idx]]+= [idx]

	#####
	d= 0
	ref_pop_dict= {}
	ssamp_dict= {}
	count_data= {}

	count_sims= list(sim_dict.keys())
	print(count_sims)
	for sim,sim_idx in sim_dict.items():

		props, pop_list, size_file, pop_file= sim_countPrep(count_array, info_array, sim, sim_idx)

		if not len(props):
			continue
		pop_diffs, ref_props, pop_props, pop_refkeys, pop_size_idx_dict, pop_size_dict, pop_idx_dict= stats_condense(props, pop_list, size_file, pop_file)
		


		for pop,si_dict in pop_size_idx_dict.items():
			print(sim,pop)
			if pop not in ref_pop_dict.keys():
				ref_pop_dict[pop]= [ref_props[pop]]
			else:
				ref_pop_dict[pop].append(ref_props[pop])

			if pop not in ssamp_dict.keys():
				ssamp_dict[pop]= {}

			for si,si_idx in si_dict.items():
				if si not in ssamp_dict[pop].keys():
					ssamp_dict[pop][si]= []

				for idx in si_idx:
				    #local_array= count_array[pop_idx_dict[pop][idx],:]
				    ssamp_dict[pop][si].append(pop_props[pop][idx])#.reshape(row,col))
	counts_dict= {}
	stats_dict= {}

	for pop in ssamp_dict.keys():
		counts_dict[pop]= {'sizes':[]}
		stats_dict[pop]= {}
		for size in sorted(ssamp_dict[pop].keys()):
			if size > si_max:
				continue
			counts_dict[pop]['sizes'].append(size)
			set1= ssamp_dict[pop][size]

			dists_self= set_SSD(set1,set1,same= True)
			dists_ref= set_SSD(set1,ref_pop_dict[pop],same= False)
			dists_ref= np.array(dists_ref).reshape(len(set1),len(ref_pop_dict[pop]))
			
			dists_ref= np.mean(np.array(dists_ref),axis= 1)

			stats_dict[pop][size]= {
				'self': dists_self,
				'ref': dists_ref
				}
		print("##### reff pops")
		print(ref_pop_dict[pop])
		ref_pop_dict[pop]= set_SSD(ref_pop_dict[pop],ref_pop_dict[pop],same= True)
		print(ref_pop_dict[pop])

	return stats_dict, counts_dict, ref_pop_dict


