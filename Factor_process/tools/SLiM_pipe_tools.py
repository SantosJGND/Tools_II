import sys
import argparse
import numpy as np

from itertools import product
import itertools as it

import tempfile
import os
import gzip 
import subprocess

from datetime import datetime




def write_popIDs(sampleSizes,popSuff= 'pop',indSuff= 'i',
                file_dir= '/mnt/d/GitHub/fine-scale-mutation-spectrum-master/slim_pipe/mutation_counter/data/'):
    '''create identifiers file ind/pops'''
    
    popNames= [popSuff + str(x) for x in range(len(sampleSizes))]
    indIDs= [indSuff + str(x) for x in range(sum(sampleSizes))]
    
    popIDs= np.repeat(np.array(popNames),sampleSizes)
    data= np.array([
        indIDs,
        popIDs
    ]).T
    
    filename= file_dir + '/ind_assignments.txt'
    
    with open(filename,'w') as f:
        vector= ['\t'.join(x) for x in data]
        vector= '\n'.join(vector)
        f.write(vector)



def read_chrom_sizes(assembly,size_dir= 'chrom_sizes/',reject= ['M','Un','X','Y','_']):
    '''
    read chromosome size file. Store and return as {str(chrom): int(chrom_len)} dict.
    '''
    filename= size_dir + assembly + '.chrom.sizes'
    with open(filename,'r') as f:
        lines= f.readlines()
        lines= [line.strip().split() for line in lines]
        
        Sizes= {x[0].strip('chr'): int(x[1]) for x in lines if sum([int(y in x[0]) for y in reject]) == 0}
    
    return Sizes



def region_samplev2(L, chrom_sizes, N, fasta_file= ''):
    '''
    prepare sequence dictionary: {chrom:list(start pos)}
    provide to fasta_Rextract.
    
    ##
    '''
    
    chroms= list(chrom_sizes.keys())
    sizes= [int(chrom_sizes[x]) for x in chroms]
    sizes= np.array(sizes) / sum(sizes)
    
    choices= np.random.choice(chroms,size= N,p= sizes,replace= True)
    
    seqs_store= {
        z: len([x for x in range(N) if choices[x]==z]) for z in list(set(choices))
    }
    
    seqs= fasta_RextractUnif(fasta_file, seqs_store, L= L)
    
    return seqs



def fasta_RextractUnif(fasta,seq_store,L= 10000):
    ''' 
    Extract regions from fasta file.
    Takes dictionary {chrom: list(start points)};
    extracts sequences after reading each chromosome.
    '''
    refseqs= {x:{} for x in seq_store.keys()}
    
    Outfiles= {
        x: tempfile.TemporaryFile() for x in seq_store.keys()
    }
    
    d= 0
    
    with gzip.open(fasta,'rb') as fp:
        for line in fp:
            #line= line.decode()
            if line[0:1] == b'>':
                head= line.decode()
                head= head.strip()
                if d != 0:
                    d=0

                head= head[1:].split('\t')
                head= head[0].strip('chr')

                if head in seq_store.keys():
                    print('opening fasta chr: {}'.format(head))
                    d= head
                    outfile= Outfiles[d]
                    continue
            
            if d != 0:
                processed_line= line.upper().strip(b'\n')
                outfile.write(processed_line)
    
    for chrom in Outfiles.keys():
        f= Outfiles[chrom]
        f.seek(os.SEEK_SET)
        result= f.read().decode()
        f.close()
        
        chrom_seqs= return_seqs(result,size= seq_store[chrom],L= L)
        
        refseqs[chrom]= chrom_seqs
    
    return refseqs



def return_seqs(seq,size= 10,L= 1000,keep= ['A','T','G','C'],threshold= 0.01):
    '''returns N=size segments of length L, unique.()=keep'''
    d= 0
    
    seqL= len(seq)
    seq_dict= {}

    lt= threshold * L
    
    while d < size:
        pos= np.random.randint(low= 0, high= seqL - L,size= 1)[0]
        given= seq[pos:(pos + L)].upper()
        
        scrag= [x for x in range(len(given)) if given[x] not in keep]
        
        if len(scrag) < lt:
            if scrag:
                replace= np.random.choice(keep,len(scrag))
                given= list(given)
                for idx in range(len(scrag)):
                    given[scrag[idx]] = replace[idx]
                given= ''.join(given)

            seq_dict[pos]= given
            
            d += 1
    
    return seq_dict



def write_fastaEx(fasta,chrom= '1',start= 0, ID= 'SIM', fasta_dir= ''):
    ''' write extracted sequence as fasta file'''
    filename= ''.join([fasta_dir,'/chr',chrom,'_',ID,'.fa'])
    
    header= ' '.join(['>'+ID,chrom,str(start),str(start + len(fasta))]) + '\n'
    
    with open(filename,'w') as fp:
        fp.write(header)
        fp.write(fasta)
    
    return filename




def process_recipe(recipe,constant_dict, SIMname, remove_dirs= False):
    '''add constant definitions to a SLiM recipe'''
    
    new_recipe= recipe.split('/')
    new_recipe[-1]= SIMname + '_' + new_recipe[-1]
    new_recipe= '/'.join(new_recipe)
    
    with open(recipe,'r') as f:
        lines= f.readlines()
    
    init_idx= [x for x in range(len(lines)) if 'initialize()' in lines[x]][0]
    
    defConst= '\tdefineConstant({},{});\n'
    
    for v,g in constant_dict.items():
        if v != "other":
            if isinstance(g,str):
                if remove_dirs:
                    g= g.split('/')[-1]
                newline= defConst.format('"{}"'.format(v),'"{}"'.format(g))
            else:
                newline= defConst.format('"{}"'.format(v),str(g))
            
            lines.insert(init_idx+1,newline)
    
    with open(new_recipe,'w') as f:
        f.write(''.join(lines))
    
    return new_recipe


def process_recipe_other(recipe,constant_dict, SIMname):
    '''add constant definitions to a SLiM recipe'''
    
    with open(recipe,'r') as f:
        lines= f.readlines()
    
    other= constant_dict["other"]
    other_idx= {
        z: [x +1 for x in range(len(lines)) if z in lines[x]][0] for z in other.keys()
    }
    
    for v in other_idx.keys():
        lines.insert(other_idx[v],other[v])

    with open(recipe,'w') as f:
        f.write(''.join(lines))
    
    return recipe





def SLiM_dispenserv1(sim_store, sim_recipe, cookID= 'ID', slim_dir= './', batch_name= '',
                    ID= 'v1',L= 10000, logSims= 'sims.log', mutlog= 'toMut.log'):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    '''
    nf= len(sim_store)
    for SIMname in sim_store.keys():
        
        command_line_constants= sim_store[SIMname]
        
        ### generate modified slim recipe
        new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname)

        if "other" in command_line_constants:
            new_recipe= process_recipe_other(new_recipe,command_line_constants,SIMname)
        
        seed= np.random.randint(0,high= nf,size= 1)[0]
        ### Launch SLiM through shell.
        slim_soft= slim_dir + 'slim' 
        
        command_units= [slim_soft, '-m', '-s', str(seed), new_recipe]
        command_units= ' '.join(command_units)
        print(command_units)
        os.system(command_units)

        os.system('gzip {}'.format(sim_store[SIMname]["vcf_file"]))
        os.system('gzip {}'.format(sim_store[SIMname]["fasta_file"]))
        os.system('rm {}'.format(new_recipe))

        now = datetime.now()
        tnow= now.strftime("%d/%m/%Y %H:%M:%S")
        
        constant_str= ';'.join(['='.join([v,str(g)]) for v,g in command_line_constants.items()])

        with open(logSims,'a') as fp:
            INFO= [SIMname,tnow,sim_recipe,cookID,'L='+str(L),constant_str]
            fp.write('\t'.join(INFO) + '\n')
        
        with open(mutlog,'a') as fp:
            fp.write(SIMname + '\n')



def SLiM_dispenserv2(sim_store, cookID= 'ID', slim_dir= './', batch_name= '',
                    ID= 'v1',L= 10000, logSims= 'sims.log', mutlog= 'toMut.log'):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    - v2 finds recipes within sim_store instead of using same global recipe. 
    '''
    nf= len(sim_store)
    for SIMname in sim_store.keys():
        
        command_line_constants= sim_store[SIMname]
        sim_recipe= command_line_constants['recipe']

        command_line_constants.pop('recipe')
        ### generate modified slim recipe
        new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname)

        if "other" in command_line_constants:
            new_recipe= process_recipe_other(new_recipe,command_line_constants,SIMname)
        
        seed= np.random.randint(0,high= nf,size= 1)[0]
        ### Launch SLiM through shell.
        slim_soft= slim_dir + 'slim' 
        
        command_units= [slim_soft, '-m', '-s', str(seed), new_recipe]
        command_units= ' '.join(command_units)
        
        os.system(command_units)

        os.system('gzip {}'.format(sim_store[SIMname]["vcf_file"]))
        os.system('gzip {}'.format(sim_store[SIMname]["fasta_file"]))
        os.system('rm {}'.format(new_recipe))

        now = datetime.now()
        tnow= now.strftime("%d/%m/%Y %H:%M:%S")
        
        constant_str= ';'.join(['='.join([v,str(g)]) for v,g in command_line_constants.items()])

        with open(logSims,'a') as fp:
            INFO= [SIMname,tnow,sim_recipe,cookID,'L='+str(L),constant_str]
            fp.write('\t'.join(INFO) + '\n')
        
        with open(mutlog,'a') as fp:
            fp.write(SIMname + '\n')


def sbatch_launch(command,ID,batch_dir= '',modules= ['module load python/3.6.4','module load slim/3.3.1'], 
                    mem= '15GB',t= '30:00:00',nodes= 4, debug= False):
    
    header_array= [
    '''#!/bin/bash
#SBATCH -n {}
#SBATCH -t {}
#SBATCH --mem={}
module purge'''.format(nodes,t,mem),
'''#!/bin/bash
#SBATCH -t 15 -p debug -q wildfire''']
    header= header_array[int(debug)]
    lines= [header + '''\n''']
    
    for mod in modules:
        lines.append(mod + '\n')

    lines.append('\n')
    lines.append(command)

    filename= batch_dir + ID + '_batch.sh'

    with open(filename,'w') as fp:
        fp.write(''.join(lines))

    os.system('chmod +x {}'.format(filename))
    
    return filename




def SLiM_dispenserv3(sim_store, sim_recipe= '', cookID= 'ID', slim_dir= './', batch_name= '',
                    ID= 'v1',L= 10000, logSims= 'sims.log', mutlog= 'toMut.log', mem= '15GB',t= '30:00:00',nodes= 4):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    '''
    nf= len(sim_store)
    for SIMname in sim_store.keys():
        command_line_constants= sim_store[SIMname]

        sim_dir= command_line_constants["vcf_file"].split('/')[:-1]
        sim_dir= '/'.join(sim_dir) + '/'


        if not sim_recipe:
            sim_recipe= command_line_constants['recipe']
            command_line_constants.pop('recipe')
        
        ### generate modified slim recipe
        new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname)

        if "other" in command_line_constants:
            new_recipe= process_recipe_other(new_recipe,command_line_constants,SIMname)
        
        seed= np.random.randint(0,high= nf,size= 1)[0]
        ### Launch SLiM through shell.
        slim_soft= slim_dir + 'slim' 
        
        command_units= [slim_soft, '-m', '-s', str(seed), new_recipe]
        command_units= ' '.join(command_units)
        print(command_units)

        sbatch_file= sbatch_launch(command_units,SIMname,batch_dir= sim_dir, mem= mem,t= t,nodes= nodes)

        os.system('sbatch ' + sbatch_file)

        now = datetime.now()
        tnow= now.strftime("%d/%m/%Y %H:%M:%S")
        
        constant_str= ';'.join(['='.join([v,str(g)]) for v,g in command_line_constants.items()])

        with open(logSims,'a') as fp:
            INFO= [SIMname,tnow,sim_recipe,cookID,'L='+str(L),constant_str]
            fp.write('\t'.join(INFO) + '\n')
        
        with open(mutlog,'a') as fp:
            fp.write(SIMname + '\n')




def osg_template(osg_submit,ID,executable,input_files,output_files,arguments,
    cpus= 1,mem= 'GB',Nmem= 1,diskN= 1,diskS= 'GB',log_dir= 'log',queue= 1):

    '''
    Write osg_connect submit file as function.
    '''
    lines= []

    output_dict= {
        x: x.split('/')[-1] for x in output_files
    }

    lines.append('executable = ' + executable)
    lines.append('arguments = ' + ' '.join(arguments))
    lines.append('transfer_input_files = ' +  ','.join(input_files))
    lines.append('transfer_output_files = ' + ','.join(list(output_dict.values())))

    remaps= []
    for filepath,file in output_dict.items():
        remaps.append('='.join([file,filepath]))

    lines.append('transfer_output_remaps = "{}"'.format(' ; '.join(remaps)))

    lines.append('\n')
    lines.append('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-ubuntu-18.04:latest"')
    lines.append('Requirements = (HAS_MODULES =?= true) && (OSGVO_OS_STRING == "RHEL 7") && (HAS_SINGULARITY == TRUE)')
    lines.append('\n')

    for x in ['error','output','log']:
        lines.append('{} = {}/job.$(Cluster).$(Process).{}.{}'.format(x,log_dir,ID,x))

    lines.append('\n')
    lines.append('request_cpus = ' + str(cpus))
    lines.append('request_memory = {} {}'.format(str(Nmem),mem))
    lines.append('request_disk = {} {}'.format(str(diskN),diskS))

    lines.append('\n')
    lines.append('queue {}'.format(queue))

    with open(osg_submit,'w') as f:
        f.write('\n'.join(lines))


def SLiM_osg_dispenser(sim_store, sim_recipe= '', cookID= 'ID', slim_dir= './', batch_name= '',
                    ID= 'v1',L= 10000, logSims= 'sims.log', mutlog= 'toMut.log',dir_data= './data/',
                    cpus= 1,Nmem= 1,mem= 'GB',diskN= 1,diskS= 'GB',log_dir= 'log'):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    - v2 finds recipes within sim_store instead of using same global recipe. 
    '''
    nf= len(sim_store)
    for SIMname in sim_store.keys():
        
        command_line_constants= sim_store[SIMname]

        if not sim_recipe: 
            sim_recipe= command_line_constants['recipe']
            command_line_constants.pop('recipe')
        

        ### generate modified slim recipe
        new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname,remove_dirs= True)
        rec_stdlone= new_recipe.split('/')[-1]
        #rec_stdlone= '/'.join(rec_stdlone)

        if "other" in command_line_constants:
            new_recipe= process_recipe_other(new_recipe,command_line_constants,SIMname)
        
        seed= np.random.randint(0,high= nf,size= 1)[0]
        ### Launch SLiM through shell.
        slim_soft= slim_dir + 'slim' 
        
        osg_program= slim_soft
        osg_arguments= ['-m', '-s', str(seed), rec_stdlone]
        osg_files= [new_recipe,command_line_constants['fasta_file']]

        if 'mut_file' in command_line_constants.keys():
            osg_files.append(command_line_constants['mut_file'])
        
        osg_output= [command_line_constants["vcf_file"]]

        osg_submit= dir_data + SIMname + '/' + SIMname + '.submit'

        osg_template(osg_submit,SIMname,osg_program,osg_files,osg_output,osg_arguments,cpus= cpus,Nmem= Nmem,mem= mem,diskN= diskN,diskS= diskS)

        os.system('condor_submit ' + osg_submit)

        now = datetime.now()
        tnow= now.strftime("%d/%m/%Y %H:%M:%S")
        
        constant_str= ';'.join(['='.join([v,str(g)]) for v,g in command_line_constants.items()])

        with open(logSims,'a') as fp:
            INFO= [SIMname,tnow,sim_recipe,cookID,'L='+str(L),constant_str]
            fp.write('\t'.join(INFO) + '\n')
        
        with open(mutlog,'a') as fp:
            fp.write(SIMname + '\n')




