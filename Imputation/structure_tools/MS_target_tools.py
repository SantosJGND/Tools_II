
def MS_get_norm(Sequences,refs_lib,ncomps= 4,clsize= 15,Bandwidth_split= 20,
               pca_qtl= 0.2):
    '''
    Perform PCA + Mean Shift across windows. Extract Meanshift p-value vectors. Perform amova (optional).
    '''

    pca = PCA(n_components=ncomps, whiten=False,svd_solver='randomized').fit(Sequences)
    data = pca.transform(Sequences)

    params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0,cv= 3, iid= False)

    ######################################
    ####### TEST global Likelihood #######
    ######################################
    Focus_labels = [z for z in it.chain(*refs_lib.values())]

    #### Mean Shift approach
    ## from sklearn.cluster import MeanShift, estimate_bandwidth

    bandwidth = estimate_bandwidth(data, quantile= pca_qtl, n_samples=len(Focus_labels))
    if bandwidth <= 1e-3:
        bandwidth = 0.1

    ms = MeanShift(bandwidth=bandwidth, cluster_all=False, min_bin_freq=clsize)
    ms.fit(data[Focus_labels,:])
    labels = ms.labels_


    Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}
    Keep= [x for x in Tree.keys() if len(Tree[x]) > clsize]

    Tree= {x:Tree[x] for x in Keep}
    Ngps= len(Tree)

    ### Extract MScluster likelihood by sample

    dist_store= {}

    for hill in Tree.keys():
        
        grid.fit(data[Tree[hill],:])

        # use the best estimator to compute the kernel density estimate
        kde = grid.best_estimator_

        # normalize kde derived log-likelihoods, derive sample p-values
        P_dist = kde.score_samples(data[Tree[hill],:])
        Dist = kde.score_samples(data)
        P_dist= np.nan_to_num(P_dist)
        Dist= np.nan_to_num(Dist)
        
        if np.std(P_dist) == 0:
            Dist= np.array([int(Dist[x] in P_dist) for x in range(len(Dist))])
        else:
            Dist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Dist)
            Dist= np.nan_to_num(Dist)
            dist_store[hill]= Dist
    
    return Tree, dist_store,data



def kde_gen_dict(data,label_dict):
    '''
    create dictionary of group kde generators in data space.
    '''
    
    params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0,cv= 3, iid= False)

    ref_gens= {}
    ref_stats= {}

    for hill in label_dict.keys():

        grid.fit(data[label_dict[hill],:])
        # use the best estimator to compute the kernel density estimate
        kde = grid.best_estimator_
        ref_gens[hill]= kde
        
        kd_scores= kde.score_samples(data[label_dict[hill],:])
        kd_stats= [np.mean(kd_scores),np.std(kd_scores)]
        ref_stats[hill]= kd_stats
    
    return ref_gens, ref_stats

from scipy.stats import norm



def gen_class(samples,ref_generators,gen_stats= {},lb= 1e-3,out_code= -1):
    '''
    use kde generators in dictionary to score and classify samples.
    '''
    ref_keys= list(ref_generators.keys())
    score_dict= {z: g.score_samples(samples) for z,g in ref_generators.items()}
    if gen_stats:
        
        score_dict= {z: norm.cdf(g,loc= gen_stats[z][0],scale= gen_stats[z][1]) for z,g in score_dict.items()}
    #print([x.shape for x in score_dict.values()])
    score_array= [score_dict[z] for z in ref_keys]
    score_array= np.array(score_array)
    #score_array= np.exp(score_array)
    
    maxs= np.max(score_array,axis= 0)
    #print(maxs)
    maxs= maxs < lb
    
    score_sum= np.sum(score_array,axis= 0)
    score_sum[score_sum == 0]= 1
    score_array= score_array / score_sum
    
    maxl= np.argmax(score_array,axis= 0)

    maxl= np.array(ref_keys)[maxl]
    maxl[maxs]= out_code
    
    return maxl


def clustClass(ms_local,pca_obj,ref_gens,gen_stats= {},out_code= -1, 
               return_mean= True,lb= 1e-2):
    '''
    ms_local= distances by cluster.
    '''
    
    mskeys= list(ms_local.keys())

    ## 
    dist_array= [ms_local[g] for g in mskeys]
    dist_array= np.array(dist_array)
    qtl_dist= pca_obj.transform(dist_array)
    #print(qtl_dist.shape)
    ## Classify kde profiles. 
    cluster_class= gen_class(qtl_dist,ref_gens,gen_stats= gen_stats,lb= lb, 
                             out_code= out_code)
    
    
    cluster_found= {z: [x for x in range(len(cluster_class)) if cluster_class[x] == z] for z in list(set(cluster_class)) if z != out_code}

    for v,g in cluster_found.items():
        dist_foud= qtl_dist[g]
        if dist_foud.shape[0] > 1:
            dist_foud= np.mean(dist_foud,axis= 1)

        g= dist_foud    
    
    return cluster_found
