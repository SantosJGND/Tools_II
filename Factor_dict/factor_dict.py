## 

import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


def array_to_dict_hierarch(vector_list,sect_array= [],set_dict= {},n= 0):
    '''
    Convert list of factor vectors into to nested set dictionary of indices.
    Vectors must be of the same size.
    '''
    
    if n== 0:
        sect_array= list(range(len(vector_list[0])))
        set_dict= {
            0: sect_array
        }
    
    if n == len(vector_list):
        return sect_array
        ## or do something else. like, plotting only using these values.
        ## for example by storing list of nodes travelled as title or filename. 
    
    else:
        parent= vector_list[n]
    
        indices= sect_array
        #
        
        vector= [parent[x] for x in indices]
        
        subset_dict= recursively_default_dict()

        for idx in range(len(vector)):
            subset_dict[vector[idx]][idx]= 0

        subset_dict= {
            z: sorted(list(w.keys())) for z,w in subset_dict.items()
        }
        
        subset_dict= {
            z: [indices[x] for x in g] for z,g in subset_dict.items()
        }
        
        
        return {
            z: array_to_dict_hierarch(vector_list,sect_array= g,n= n+1) for z,g in subset_dict.items()
        }


#hierarchy= [factor1,factor2,..,factorN]
#batch_pop_ref_dict=  array_to_dict_hierarch(hierarchy,sect_array= [],set_dict= {},n= 0)

