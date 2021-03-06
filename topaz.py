'''
topaz.py package
Assembles miRNA, transcriptional data from garnet and gene expression data into one single network
to identify putative transcriptional regulations

Requires SAMNet to run network algorithm and GARNET to process epigenetic data

Copyright (c) 2014-2015 Sara JC Gosline

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''
__author__='Sara JC Gosline'
__email__='sgosline@mit.edu'

import pickle,re,os,sys,random
import networkx as nx
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from copy import deepcopy
from collections import defaultdict

##need to update this if samnet is installed somewhere else
#sys.path.append("../SAMNet/src")
#global samnet


from itertools import tee, izip
def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

##all levels of hte graph (from top to bottom)
all_graph_levels=['mirAbund','mirFC','mirTarget','TFexpr','motif','mRNA']

def build_samnet_network(miRVals,miRFcs,miRTargs,upRegulatedGenes,downRegulatedGenes,tfsOfInterest,upRegulatedChromRegions,downRegulatedChromRegions,do_hier=True,delim='.',cutoff=0.5,draw=False): #,use_down_mirs=True):
                         #top_mirs=15,upstream='',thresh=0.5,do_hier=True, default_node_cap=0.0,use_clip=True,use_cds_clip=False,dr_tfs={},ko_mrnas=[],wt_mrnas=[],ppi_file='',ppi_mapping_file='',hmarks=['H3K4me3','H3K27ac'],largeClusters=True,donorm=False):
    '''
    Assembles miRNAs and TFs and histone data into network to enable running by SAMNet
    miRVals: list of miRNAs and their abundance or FC
    miRFcs: list of fold changes, if available (can be empty)
    miRTargs: dictionary of dictionary of miRNAs, targets and scores
    upRegulatedGenes: list of up-regulated genes
    downRegulatedGenes: list of down-regulated genes.
    tfsOfInterest: list of transcription factors of interest
    upRegulatedChromRegions: dictionary of TF-mRNA networks. Keys indiciate histone mark
    downRegulatedChromRegions: dictionary of TF-mRNA networks. Keys indiciate histone mark
    do_hier: add in hierarchical edge capacities
    delim: delimiter used in GARNET To collapse TF names
    cutoff: score below which we don't include TF-DNA interactions
    draw: try to draw network
    '''
    ##first start with all histone marks, and add in chrom-region-mRNA level of network
    hmarks=upRegulatedChromRegions.keys()

    #distinguish TFs by whether or not they are up-regulated of down-regulated
    upperCaseTfs={}
    for k in tfsOfInterest.keys():
        upperCaseTfs[k.upper()]=tfsOfInterest[k]

    up_tfs=[a for a in upperCaseTfs.keys() if upperCaseTfs[a]>0]
    down_tfs=[a for a in upperCaseTfs.keys() if upperCaseTfs[a]<0]
    print 'Have '+str(len(up_tfs))+' up-regulated TFs and '+str(len(down_tfs))+' down-regulated'
    

    #distinguish miRNAs similarly
    if len(miRFcs)==0:
        up_mirs=[a for a in miRVals.keys() if miRVals[a]>0]
        down_mirs=[a for a in miRVals.keys() if miRVals[a]<0]
    else:
        up_mirs=[a for a in miRFcs.keys() if miRFcs[a]>0]
        down_mirs=[a for a in miRFcs.keys() if miRFcs[a]<0]
        
    print 'Have '+str(len(up_mirs))+' up-regulated miRs and '+str(len(down_mirs))+' down-regulated'
    
    #assemble weighted graph
    graph=nx.DiGraph()
    
    ##create reverse targets in case we have fold change AND abundance values
    reverseTargs=defaultdict(dict)
    for mir in miRTargs.keys():
        for key,val in miRTargs[mir].iteritems():
            reverseTargs[key][mir]=val
    mir_to_mirfam=defaultdict(list)
    mirfam_to_targs=defaultdict(dict)
    for tf in reverseTargs.keys():
        scoredict=defaultdict(list)
        for mir,val in reverseTargs[tf].iteritems():
            scoredict[val].append(mir)
        
        for score,mirs in scoredict.iteritems():
            if len(mirs)==1: #create dummy node to account for fold change
                joinednode=mirs[0]+'.'+mirs[0]
            else:
                joinednode='.'.join(mirs)
            mirfam_to_targs[joinednode][tf]=score
            for mir in mirs:
                mir_to_mirfam[mir].append(joinednode)
                
    #now start with miRNas
    up_tfs_in_play=set()
    down_tfs_in_play=set()
    for m in miRVals.keys():
        if m not in miRTargs.keys():
            print 'miRNA %s has no target information, skipping'%(m)
            continue

        #print 'Evaluating %d targets for miR %s:'%(len(targs),m)
        #print ','.join(targs)
        ##allow up-regulation and down-regulation to be explained by miRNAs in either direction
        graph.add_edge('Up',m,weight=1.0)###remove this node at end, just a dummy for sp computation
        graph.add_edge('Down',m,weight=1.0)###remove this node at end, just a dummy for sp computation

        if len(miRFcs)==0:
            #print 'Connecting miRNAs directly to target TFs'
            targs=miRTargs[m].keys()
            for targ in targs:
                w=np.abs(float(miRTargs[m][targ]))
                if m in up_mirs and targ.upper() in down_tfs and not np.isnan(w):
                    graph.add_edge(m,targ.upper(),weight=np.abs(w))
                    up_tfs_in_play.add(targ.upper())
                elif m in down_mirs and targ.upper() in up_tfs and not np.isnan(w):
                    graph.add_edge(m,targ.upper(),weight=np.abs(w))
                    down_tfs_in_play.add(targ.upper())
        else:
            try:    
                fcscore=np.abs(miRFcs[m])
            except:
                print 'No Fold change value for %s, removing from network'%(m)
                continue
            #print 'Treating FC and Abundance as separate values through miRNA score families'
            for mirfam in mir_to_mirfam[m]:
                for targ,val in mirfam_to_targs[mirfam].iteritems():
                    targscore=np.abs(float(val)) #get target score
                    if np.isnan(targscore):
                        continue
                    if m in up_mirs and targ.upper() in down_tfs:
                        graph.add_edge(m,mirfam,weight=fcscore)
                        graph.add_edge(mirfam,targ.upper(),weight=targscore)
                        up_tfs_in_play.add(targ.upper())
                    elif m in down_mirs and targ.upper() in up_tfs:
                        graph.add_edge(m,mirfam,weight=fcscore)
                        graph.add_edge(mirfam,targ.upper(),weight=targscore)
                        down_tfs_in_play.add(targ.upper())

    ##if miRNA values are not changing, the network will be empty at this point. So instead 
    print 'Graph has '+str(len(graph.nodes()))+' nodes and '+str(len(graph.edges()))+' edges from '+str(len(miRTargs))+' possible miRNAs to '+str(len(up_tfs_in_play))+' possible up-regulated TFs and '+str(len(down_tfs_in_play))+' downregulated TFs'
    
    all_tfs_in_play=[a for a in up_tfs_in_play]+[a for a in down_tfs_in_play]
    ##then create TF-cond.Tf.histone interactions
    ###first start with Up######################################
    
    cond='Up'
    tfs_with_binding=set()
    mrnas=set()
    for hist in upRegulatedChromRegions.keys():
        #for every histone create individual node for each matrix
        chromGraph=upRegulatedChromRegions[hist]
        for matrix,targs in chromGraph.adjacency_iter():
            if len(targs)==0:
                continue
            #figure out which TFs target this matrix
            #print matrix.split(delim)
            tfs=[t.upper() for t in matrix.split(delim) if t.upper() in all_tfs_in_play] #allow only up-regulated TFs
            #print tfs
            ##add edge from every TF to the matrix
            has_binding=False
            for a in tfs:
                weight=np.abs(upperCaseTfs[a])##take absolute value of fold change, relying on miRNAs to filter
#                print 'Weight from '+a+' to '+matrix+' is :'+str(weight)
                if weight>0.0:
                    graph.add_edge(a,hist+'.'+cond+'.'+matrix,weight=weight)
                    has_binding=True
                tfs_with_binding.add(a)
            if not has_binding:
                continue
            #now iterate through every target and add edge from matrix to target
            for targ,eattr in targs.items():
                weight=eattr['weight']
                targ=re.sub("mrna",'',targ)##in case we added it
#                print targ,str(weight)
                if weight>cutoff and targ in upRegulatedGenes:
                    graph.add_edge(hist+'.'+cond+'.'+matrix,targ+'mrna',{'weight':weight})
                    mrnas.add(targ)
                   # print 'adding up target '+targ
    #        print 'Graph has '+str(len(graph.nodes()))+' and '+str(len(graph.edges()))+' after '+hist+' Up regions and mrnas added'
    print 'TFs in up-reg regions: '+','.join([a for a in tfs_with_binding])

    cond='Down'
    for hist in downRegulatedChromRegions.keys():
        #for every histone create individual node for each matrix
        chromGraph=downRegulatedChromRegions[hist]
        for matrix,targs in chromGraph.adjacency_iter():
            if len(targs)==0:
                continue
            #figure out which TFs target this matrix
#            if len(down_mirs)==0:##in this case we let up-regulated TFs down-regulate?
            tfs=[t.upper() for t in matrix.split(delim) if t.upper() in all_tfs_in_play]
#            else:
#                tfs=[t.upper() for t in matrix.split(delim) if t.upper() in down_tfs_in_play]
            #print tfs
            ##add edge from every TF to the matrix
            has_binding=False
            for a in tfs:
                weight=np.abs(upperCaseTfs[a])##take absolute value of fold change, relying on miRNAs to filter
                if weight>0.0:
                    graph.add_edge(a,hist+'.'+cond+'.'+matrix,weight=weight)
                    has_binding=True
                tfs_with_binding.add(a)
            #now iterate through every target and add edge from matrix to target
            if not has_binding:
                continue
            for targ,eattr in targs.items():
                weight=eattr['weight']
                targ=re.sub("mrna",'',targ)##in case we added it
                 #               print targ,str(weight)
                if weight>cutoff and targ in downRegulatedGenes:
                    graph.add_edge(hist+'.'+cond+'.'+matrix,targ+'mrna',{'weight':weight})
                    mrnas.add(targ)
                   # print 'adding up target '+targ
    #        print 'Graph has '+str(len(graph.nodes()))+' and '+str(len(graph.edges()))+' after '+hist+' Down regions and mrnas added'
    print 'TFs in all regions: '+','.join([a for a in tfs_with_binding])

    tfs_in_play=[a for a in up_tfs_in_play]+[a for a in down_tfs_in_play]
    print str(len(mrnas))+' mRNAs in network'
    print 'Only '+str(len(tfs_with_binding))+' tfs have binding out of '+str(len(tfs_in_play))+' mRNA above activity threshold, removing extra nodes'

    for n in tfs_in_play:
        if n not in tfs_with_binding:
            graph.remove_node(n)
            
    ##now collect node capacities -- this will be blank if we don't build custom capacities
    node_caps=defaultdict(dict)
    if do_hier:
        print 'Constructing node-specific capacties'
        ##Hierarchical capacities more fully mimic the flow model - miRNA-TF edges can handle more flow then TF-mRNA capacities...
        for node in graph.nodes():
            for treat in ['Up','Down']:
                if node!=treat and node!=treat+'_sink':
                    spdist=1.0
                    try:
                        spdist=nx.shortest_path_length(graph,treat,node,None)
                    except nx.exception.NetworkXNoPath:
                        continue
                    node_cap=np.power(10.0,float(spdist)*-1.0)
                    #print treat,node,str(spdist),str(node_cap)
                    node_caps[treat][node]=node_cap

    upsuc=graph.successors('Up')
    downsuc=graph.successors('Down')

    if 'Up' in graph.nodes():
        graph.remove_node('Up')
    if 'Down' in graph.nodes():
        graph.remove_node('Down')

    for n in graph.nodes():
        if graph.degree(n)==0:
            #print 'removing '+n
            graph.remove_node(n)
        
    #print 'Source nodes for up commodity '+','.join(upsuc)
    #print 'Source nodes for down commodity '+','.join(downsuc)
            
    #pickle.dump(graph,file=open('graph.pkl','w'))
    if draw:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        try:
            from networkx import graphviz_layout
            #layout=graphviz_layout
        except ImportError:
            raise('This example needs Graphvis and either Pygraphvis or Pydot')
            #layout=nx.spring_layout
    
        plt.figure(figsize=(8,8))
        nodetypes={'mir':'g','tf':'y','hist':'b','mrna':'r'}#]
        nodecols=[]
        for node in graph:
            ntype='hist'
            if 'mrna' in node:
                ntype='mrna'
            if node in tfs_in_play:
                ntype='tf'
            elif node in miRVals.keys():
                ntype='mir'
            nodecols.append(nodetypes[ntype])#nodetypes.index(ntype))
        print 'Assigned '+str(len(nodecols))+' node colors'
        pos=nx.spring_layout(graph)#graphviz_layout(graph,prog='neato')#,prog='dot',args='')
        nx.draw(graph,pos,node_color=nodecols,alpha=0.4,edge_color='grey',with_labels=False,node_size=25)

    #not sure what this does
        xmax=1.02*max(xx for xx,yy in pos.values())
        ymax=1.02*max(yy for xx,yy in pos.values())
        plt.xlim(-.02,xmax)
        plt.ylim(-.02,ymax)
        try:
            plt.legend()
        except Error:
            print 'error'
        plt.savefig('mir_tf_histone_graph.png')

        
    return graph,node_caps


def perturbNetwork(network,mirlist=[],ishier=True):
    '''
    Assume hierarchical network for now. Provide list of input weights to initiate hierarchy
    '''

    newnetwork=defaultdict(dict)
    oldnetwork=nx.to_dict_of_dicts(network)
    ##collect, for each level, all members of subsequent level
    tfs,motifs=set(),set()

    #first, collect all weights
    mirweights=[]
    for m in mirlist:
        if m not in oldnetwork.keys():
            #print 'No %s in network'%(m)
            continue
        allweights=[oldnetwork[m][d]['weight'] for d in oldnetwork[m].keys()]
        tfs.update(oldnetwork[m].keys())
        mirweights.extend(allweights)
     #   print 'Have %d edge weights for miRNA %s'%(len(allweights),m)
    np.random.shuffle(mirweights)

    print mirweights
    
    count=0
    for m in mirlist:
        if m not in oldnetwork.keys():
            #print 'No %s in network'%(m)
            continue
        for tf in oldnetwork[m].keys():
            newnetwork[m][tf]={'weight':mirweights[count]}
            count+=1
   # print 'Have %d nodes in new network'%(len(newnetwork.keys()))
            
    #now onto tfs
    tfweights=[]
    for t in tfs:
        allweights=[oldnetwork[t][d]['weight'] for d in oldnetwork[t].keys()]
        motifs.update(oldnetwork[t].keys())
        tfweights.extend(allweights)
    #    print 'Have %d edge weights for TF %s'%(len(allweights),t)
        
    np.random.shuffle(tfweights)
    count=0
    for t in tfs:
        for mot in oldnetwork[t].keys():
            newnetwork[t][mot]={'weight':tfweights[count]}
            count+=1

    #print 'Have %d nodes in new network'%(len(newnetwork.keys()))

    #now motifs, then we're done
    motweights=[]
    for mo in motifs:
        allweights=[oldnetwork[mo][d]['weight'] for d in oldnetwork[mo].keys()]
        motweights.extend(allweights)
      #  print 'Have %d edge weights for Motif %s'%(len(allweights),mo)
        
    np.random.shuffle(motweights)
    count=0
    for mo in motifs:
        for mr in oldnetwork[mo].keys():
            newnetwork[mo][mr]={'weight':motweights[count]}
            count+=1

    print 'Have %d nodes in new network'%(len(newnetwork.keys()))
    
    return nx.DiGraph(newnetwork)


def add_flow_dict(origd,newd):
    #adds flow to new dictionary
    for p,val in newd.iteritems():
        try:
            origd[p]+=val
        except:
            origd[p]=val
    return origd


def make_cyto28_files(mrna_file,mirna_file,sif_file,outdir,mirclust=False):
    '''
    This method will make cytoscape 2.8 node attribute files out of mrna and mirna expression inputs
    As well as collapse the histones out of the sif output
    '''
    allmirs=set()
    allmrnas=set()
    newfile=open(outdir+'_expressionInputs.noa','w')
    newfile.write('ExpressionValues\n')
    for row in open(mrna_file,'rU').readlines():
        arr=row.split('\t')
        newfile.write(arr[0]+' = '+arr[1])
        allmrnas.add(arr[0])
    for row in open(mirna_file,'rU').readlines():
        arr=row.split('\t')
        newfile.write(arr[0]+' = '+arr[1])
        allmirs.add(arr[0])
        
    newfile.close()
    #now process siffile
    newsif=open(outdir+'_collapsedNetwork.sif','w')
    oldsif=open(sif_file,'rU').readlines()
    #first get 'tfs' out of nodetypes
    nodetfile=open(re.sub('_mcfs.sif','_node_type.noa',sif_file),'rU')
    nodetypes={}
    for row in nodetfile.readlines()[1:]:
        gene,type = row.strip().split(' = ')
        nodetypes[gene]=type

    #now we can peruse and first map tfs to genes
    mir_clust_maps=defaultdict(set)
 #   for row in oldsif:
 #       p1,ia,p2=row.strip().split()
 #       if p1=='S1' or p2=='T1':
 #           continue
 #       if p1 in allmirs and mirclust:
 #           mir_clust_maps[p2].add(p1)

 
    tf_gene_maps=defaultdict(set)
    for row in oldsif:
        p1,ia,p2=row.strip().split()
        if p1=='S1' or p2=='T1':
            continue
        if 'transcriptionfactor' in nodetypes[p2] and not mirclust: #then it's either a TF or mirCluster
                tf_gene_maps[p2].add(p1)
        elif mirclust:
            if 'protein' in nodetypes[p1] and 'transcriptionfactor' in nodetypes[p2]:
                tf_gene_maps[p2].add(p1)
            elif p1 in allmirs and len(p2.split('.'))>1:
                tf_gene_maps[p2].add(p1)
    #print tf_gene_maps
    #now repeat file scan
    for row in oldsif:
        p1,ia,p2=row.strip().split()
        if p1=='S1' or p2=='T1':# or 'transcriptionfactor' in nodetypes[p2]:
            continue
        ##now update the mappings
        if p1 in tf_gene_maps.keys():
           p1=[tf for tf in tf_gene_maps[p1]]
        else:
            p1=[p1]
        if p2 in tf_gene_maps.keys():
            continue
        for p11 in p1:
            newsif.write(p11+' '+ia+' '+p2+'\n')

    newsif.close()
            

    
def runSamNetWithMirs(network,mrna_weights,mir_weights,gamma,samnet_path,outname,conditions=['Up','Down'],leaveOut='',node_caps={},debug=False,sinkGamma=False,stepGamma=0):
    '''
    runs samnet with mirs in network - by this point mRNA nodes should have 'mrna' appended to their name to prevent direct edges...
    '''

        ##first add samnet to path
    sys.path.insert(0,samnet_path)
    sys.path.insert(0,samnet_path+'src')
    global samnet
    import samnet  ##call run_samnet command manually
    
    ##this calls the original SAMNet code - make sure you know where it is!
    flow=0.0
    #print mir_weights
    flow,phens,prots,tfs,mrnas=samnet.run_rn(PPI_with_weights=network,indirect_weights=deepcopy(mir_weights),direct_weights={},graph_tr=nx.DiGraph(),mrna_weights=mrna_weights,output=outname,updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=debug,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
    tries=0
    
    while flow==0.0 and tries<stepGamma:
        gamma=int(gamma)+1
        tries+=1
        print '...................Flow is 0, incrementing gamma by 1 to %d'%(gamma)
        flow,phens,prots,tfs,mrnas=samnet.run_rn(PPI_with_weights=network,indirect_weights=deepcopy(mir_weights),direct_weights={},graph_tr=nx.DiGraph(),mrna_weights=mrna_weights,output=outname,updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=debug,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
        
    ##phens is another word for miRNAs - keep track of everythign
    total_phens=phens
    total_tfs=tfs
    total_mrnas=mrnas
    total_prots=prots

    res_tfs=[t for t in total_prots.keys() if t not in total_tfs.keys()]
    if len(res_tfs)==0: #this occurs in cases where miRNA targets are phenotypes
        res_tfs=[t for t in total_phens.keys() if t not in total_tfs.keys()]

    print 'Res tfs '+','.join(res_tfs)


    ##collect list of output file prefixes to combine...
    outlist=[outname+'multiComm']
    mul_c='multiComm'
    ##This code will leave out single or double TFs at once
    if leaveOut.lower()=='tf':
        for tf in res_tfs:
            print 'running with '+tf+' removed'
            newnetwork=deepcopy(network)
            newnetwork.remove_node(tf)
            newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr=nx.DiGraph(),mrna_weights=mrna_weights,output=outname+'_'+tf+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
##now add flow if we observed a knockout!
            total_prots=add_flow_dict(total_prots,newprots)
            total_phens=add_flow_dict(total_phens,newphens)
            total_tfs=add_flow_dict(total_tfs,newtfs)
            total_mrnas=add_flow_dict(total_mrnas,newmrnas)
            outlist.append(outname+'_'+tf+'_REMOVED'+mul_c+'_edge_commodity.eda')
    elif leaveOut=='mir':#code can also remove miRNAs from network, haven't really explored this
        for mir in phens.keys():
            newnetwork=deepcopy(network)
            newnetwork.remove_node(mir)
            newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr=nx.DiGraph(),mrna_weights=mrna_weights,output=outname+'_'+mir+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
            #                outlist.append(outname+'_'+mir+'_REMOVED'+mul_c+'_ppi_attributes.eda')
            total_prots=add_flow_dict(total_prots,newprots)
            total_phens=add_flow_dict(total_phens,newphens)
            total_tfs=add_flow_dict(total_tfs,newtfs)
            total_mrnas=add_flow_dict(total_mrnas,newmrnas)
            outlist.append(outname+'_'+mir+'_REMOVED'+mul_c+'_edge_commodity.eda')
            print newphens

    
    #now combine outlist files into a single unified sif file
    if len(outlist)>1:
        combined=samnet.combine_single_flows_to_make_multi(outlist[1:],outlist[0],collapse_edges=True,ismcf=True)

    return total_phens,prots,total_tfs,total_mrnas


def randomizeGraphAndRun(num_iters,miRs,mirFC,mirTargs,mirWeights,mrna_weights,upreg_genes,downreg_genes,tfs,up_chroms,down_chroms,cutoff,gamma,addpath,out,conditions,orig_out,levels=all_graph_levels,countFlow=False):
    '''
    this is the primary randomization code for the Topaz framework
    num_iters: provides the number of times to randomize
    miRs: original set of miRNAs
    mirFC: mir fold change if abundances are used. can be empty
    mirTargs: original miRNA target dictionary
    mirWeights: weights on miRNAs by commodity
    mrna_weights: weights on mRNAs by commodity
    upreg_genes,downreg_genes: sets of genes for each commodity
    tfs: dictionary of TFs and fold change
    up_chroms,down_chroms: networks for both
    cutoff: histone cutoff
    gamma: gamma value 
    addpath: path to samnet
    out:
    conditions:
    orig_out:
    levels: which data should be permuted? provide a comma-delimited list
    countFlow: if this is true, only those nodes with flow values greater than original
    will be included
    '''
        ##now we can perform randomization if necessary by re-creating entire network
    print '----------------Performing '+str(num_iters)+' randomizations-----------------------'
    print 'Randomizing the following levels: '+','.join(levels)
    outlist=[]
    mircounts,protcounts,tfcounts,mrnacounts=defaultdict(int),defaultdict(int),defaultdict(int),defaultdict(int)
    mirflow,protflow,tfflow,mrnaflow=defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list)
    
    llevels=[a.lower() for a in levels]
    allTargs=tfs.keys()#set() #keep names of TFs the same
    
    #[allTargs.update(mirTargs[m].keys()) for m in mirTargs.keys()]
    print 'Have %d miRNA targets to choose from'%(len(allTargs))
    for i in range(num_iters):
        #shuffle miRNA abundances
        if 'mirabund' in llevels or ('mirfc' in llevels and len(mirFC)==0):
            print '---Shuffling miRNA abundances for iteration %d....'%(i)
            mirnames=random.sample(mirTargs.keys(),len(miRs))
            #print 'new miRNAs: '+','.join(mirnames[0:10])
            newabs=[np.abs(a) for a in miRs.values()]
            #print 'new Values: '+','.join([str(a) for a in newabs[0:10]])
#            np.random.shuffle(newabs)
            randabund=dict(zip(mirnames,newabs))
            randMirWeights={'Up':randabund,
                'Down':randabund}
            
        else:
            randabund=miRs
            randMirWeights=mirWeights
        
        #shuffle TF expression?
        if 'tfexpr' in llevels:
            print '---Shuffling TF expression for iteration %d....'%(i)
            newvals=tfs.values()
            np.random.shuffle(newvals)
            randtfs=dict(zip(allTargs,newvals))
        else:
            randtfs=tfs
            
        upTargs=[a for a in randtfs.keys() if randtfs[a]>0]
        downTargs=[a for a in randtfs.keys() if randtfs[a]<0]

        if 'mirfc' in llevels and len(mirFC)!=0:
            print '---Shuffling miRNA fold changes for iteration %d....'%(i)
            mirnames=random.sample(randabund.keys(),len(mirFC))
            #print 'new miRNAs: '+','.join(mirnames[0:10])
            newvals=[np.abs(a) for a in mirFC.values()]
            #print 'new Values: '+','.join([str(a) for a in newabs[0:10]])
#            np.random.shuffle(newabs)
            randFC=dict(zip(mirnames,newvals))
            
        else:
            randFC=mirFC

            
        #shuffle miRNA targets?
        if 'mirtarget' in llevels:
            print '---Shuffling miRNA degree and miRNA targets for iteration %d....'%(i)
            #print mirWeights
            randomMirTargs={}
            targlist=deepcopy(mirTargs.values())
            allmirs=[a for a in mirTargs.keys()]#get all mir indicies
            for m in allmirs:
                newvals=targlist.pop(random.choice(range(len(targlist))))
                sampnum=min(len(allTargs),len(newvals))
                randomMirTargs[m]=dict(zip(random.sample(allTargs,sampnum),newvals.values()[0:sampnum]))
        else:
            randomMirTargs=mirTargs
        
        #shuffle mRNA expression changes?
        if 'mrna' in llevels:
            print '---Shuffling mRNA weights for iteration %d...'%(i)
            new_mrna_weights={}
            for k in mrna_weights.keys():
                allvals=mrna_weights[k].values()
                ##keep up-regulated targets up-regulated, and down-regulated targets down-regulated
                if k in ['up','UP','Up']:
                    upreg_genes=random.sample(upTargs,len(allvals))
                    new_mrna_weights[k]=dict(zip([a+'mrna' for a in upreg_genes],allvals))
                elif k in ['down','DOWN','Down']:
                    downreg_genes=random.sample(downTargs,len(allvals))
                    new_mrna_weights[k]=dict(zip([a+'mrna' for a in downreg_genes],allvals))
                else:
                    other_genes=random.sample(allTargs,len(allvals))
                    new_mrna_weights[k]=dict(zip([a+'mrna' for a in other_genes],allvals))
                    upreg_genes=other_genes
                    downreg_genes=other_genes
        else:
            new_mrna_weights=mrna_weights
            upreg_genes=upTargs
            downreg_genes=downTargs
            
        #shuffle histone-motif events?
        if 'motif' in llevels:
            print '---Resampling motif-DNA edges for iteration %d...'%(i)
            #first colle tall possible TFs in motifs
            alltfs=set()
#            allargs=upreg_genes+downreg_genes
            for k in up_chroms.keys():
                for h,targs in up_chroms[k].adjacency_iter():
                    if len(targs)>0:
                        alltfs.update(h.split('.'))
            for k in down_chroms.keys():
                for h,targs in down_chroms[k].adjacency_iter():
                    if len(targs)>0:
                        alltfs.update(h.split('.'))
            print 'Have %d TFs in TF matrices to shuffle'%(len(alltfs))
            rand_up_chroms,rand_down_chroms={},{}
            for mark in up_chroms.keys(): #for each mark, create random graph while preserving out-degree?
                newgraph=nx.DiGraph()##add new graph
                for matrix,targs in up_chroms[mark].adjacency_iter():
                    if len(targs)==0:
                        continue
                    #WHAT SHOULD I DO HERE???
#                    newmatname='.'.join(random.sample(alltfs,len(matrix.split('.'))))
                    newmatname=matrix
                    #print matrix,newmatname
                    newtargs=random.sample(allTargs,min(len(allTargs),len(targs)))
                    tcount=0
                    for targ,eattr in targs.items():
                        if tcount<len(newtargs):
                            newgraph.add_edge(newmatname,newtargs[tcount],weight=eattr['weight'])
                        tcount+=1
                rand_up_chroms[mark]=newgraph
            
            for mark in down_chroms.keys(): #for each mark, create random graph while preserving out-degree?
                newgraph=nx.DiGraph()##add new graph
                for matrix,targs in down_chroms[mark].adjacency_iter():
                    if len(targs)==0:
                        continue
                    #newmatname='.'.join(random.sample(alltfs,len(matrix.split('.'))))
                    newmatname=matrix
                    newtargs=random.sample(allTargs,min(len(allTargs),len(targs)))
                    tcount=0
                    for targ,eattr in targs.items():
                        if tcount<len(newtargs):
                            newgraph.add_edge(newmatname,newtargs[tcount],weight=eattr['weight'])
                        tcount+=1
                rand_down_chroms[mark]=newgraph
        else:
            rand_up_chroms=up_chroms
            rand_down_chroms=down_chroms
        
        print 'Building %dth Random Network Graph...'%(i)
       
        graph,caps=build_samnet_network(randabund,randFC,randomMirTargs,upreg_genes,downreg_genes,randtfs,rand_up_chroms,rand_down_chroms,do_hier=True,cutoff=cutoff,draw=False)


#            print new_mrna_weights[k]
        print '----------------Running SAMNet: '+str(i)+'-----------------------'
        #now graph should be built, can actually run samnet!
        newout=out+'RANDOM'+str(i)
        newmirs,prots,newtfs,mrnas=runSamNetWithMirs(graph,new_mrna_weights,randMirWeights,gamma,addpath,newout,conditions=conditions,node_caps=caps,debug=False,sinkGamma=False,stepGamma=2)
        
        if os.path.exists(newout+'multiComm_edge_commodity.eda'):
            outlist.append(newout+'multiComm_edge_commodity.eda')

        for m,val in newmirs.iteritems():
            mircounts[m]+=1
            mirflow[m].append(val)
        for p,val in prots.iteritems():
            protcounts[p]+=1
            protflow[p].append(val)
        for t,val in newtfs.iteritems():
            tfcounts[t]+=1
            tfflow[t].append(val)
        for m,val in mrnas.iteritems():
            mrnacounts[m]+=1
            mrnaflow[m].append(val)
            
    combined=samnet.combine_single_flows_to_make_multi(outlist,orig_out,collapse_edges=True,ismcf=True)

    print '----------------Randomization complete, now cleaning up-----------------------'
#    outlist.reverse()
    for o in reversed(outlist):
        try:
            os.system('rm '+re.sub('multiComm_edge_commodity.eda','*',o))
        except:
            print 'File %s probably does not exist, no flow for that iteration'%(o)

    if countFlow:
        return mirflow,protflow,tfflow,mrnaflow,len(outlist)
    else:
        return mircounts,protcounts,tfcounts,mrnacounts,len(outlist)


def main():
    '''
    Parse arguments
    '''
    progdir=os.path.dirname(sys.argv[0])

    parser=OptionParser()
    ##collect mRNA expression data: tab-delimited file of differentially expressed mRNAs
    parser.add_option('--mRNAexp',dest='mRNA_file',type='string',help='Tab-delimited file containing differentially expressed genes (preferably nuclear): gene-name and log fold change. Positive values are considered up-regulated, negative values are considered down-regulated')

    ##collect miRNA values: file needs to contain miRNA family and weight
    parser.add_option('--miRNAs',dest='miRNA_file',type='string',help='Tab-delimited file containing miRNAs and weight. Positive changes will be connected to down-regulated TFs, negative changes will be connected to up-regulated TFs')

    ##NEW OPTION (2015-08-06): enable miRNA fold change to be added separately!!
    parser.add_option('--miRFoldChange',dest='mirfc_file',type='string',help="Tab-delimited file containing miRNAs and weight.  If this is included, --miRNAs file will be treated as abundances, and a separate 'level' of the network will be added to incorporate miRNA fold change in addition to miRNA abundance", default='')
    
    
    ##collect miRNA-Target data requires pickled dictionary of dictionaries
    parser.add_option('--miRNA-targets',dest='target_file',type='string',help='PKL file of dictionary of dictionaries containing miRNAs, targets and weights')
    
    ##which TFs should be considered?
    parser.add_option('--TFWeights',dest='tf_file',type='string',help='Tab-delimited file of genes to consider as regulators. It is recommended to select all genes that exhibit at least 0.5 absolute log2 fold change')
    
    ##collect TF-mRNA edges from GARNET output: data requires pickle files for each condition
    parser.add_option('--upRegulatedRegions',dest='upRegRegions_files',type='string',help='Comma-delimited list of TF-DNA networks (in PKL) for each histone or chromatin accessibility experiment unique to up-regulation')
    parser.add_option('--downRegulatedRegions',dest='downRegRegions_files',type='string',help='Comma-delimited list of TF-DNA networks (in PKL) for each histone or chromatin accessibility experiment unique to down-regulation (or same as up-regulation if you only have one set of data')
    parser.add_option('--chromatinRegionNames',dest='exp_names',type='string',help='Comma-delimited list of names of chromatin experiments, e.g. H3K4me3')
###########################
    parser.add_option('--gamma',dest='gamma',type='string',help='SAMNet gamma parameter scales number of miRNAs selected')
    parser.add_option('--doTfKO',dest='do_ko',action='store_true',default=False,help='Do TF knockdown to identify those TFs that were robust to perturation. If selected with perturb option, will perturb all edges in network.')
    parser.add_option('--perturb',dest='num_samps',type='string',default=None,help='To assess statistics of miRNA and TF selections, use this option will re-sample miRNA targets from all proteins, to assess the probability of a particular TF showing up at random')
    
    parser.add_option('--outputPrefix',dest='out',type='string',help='Prefix for output',default='samnet')
    parser.add_option('--path-to-samnet',dest='addpath',type='string',default=os.path.join(progdir,'../SAMnet'),help='To run SAMNet we require path to SAMNet module')
    parser.add_option('--draw-graph-only',dest='plot_only',default=False,action='store_true',help='Only draw graph, do not run optimization')
    parser.add_option('--histCutoff',dest='histCutoff',default='0.5',type='string',help='Minimium log fold change required for motif match to be considered')
    parser.add_option('--levelsToPerturb',dest='levels',default=','.join(all_graph_levels),help='Commma-delimited list of levels of the network, includes all by DEFAULT:%default')
    opts,args=parser.parse_args()



    print '----------------Processing arguments-----------------------'
    #get mRNAs
    mf=open(opts.mRNA_file,'rU').readlines()
    mRNAs={}
    for row in mf:
        if len(row.strip().split('\t'))!=2:
            print '--mRNAexp file needs to be tab-delimited with 2 columns, one for mRNA, one for log fold change'
            sys.exit()
        else:
            mrna,val=row.strip().split('\t')
            mRNAs[mrna]=float(val)
    print 'Read in '+str(len(mRNAs))+' mRNAs from file'

    #get miRNAs
    miRs={}
    mif=open(opts.miRNA_file,'rU').readlines()
    for row in mif:
        if len(row.strip().split('\t'))!=2:
            print '--miRNAs file needs to be tab-delimited with 2 columns, one for miR one for lfc or abundance'
            sys.exit()
        else:
            mir,val=row.strip().split('\t')
            miRs[mir.strip()]=float(val.strip())
    print 'Read in '+str(len(miRs))+' miRNA values'

    #get miRNA FCs
    mir_fcs={}
    if opts.mirfc_file!='':
        mfc=open(opts.mirfc_file,'rU').readlines()
        for row in mfc:
            if len(row.strip().split('\t'))!=2:
                print '--miRFoldChange file needs to be tab-delimited with 2 columns, one for miR one for lfc'
                sys.exit()
            else:
                mir,val=row.strip().split('\t')
                mir_fcs[mir.strip()]=float(val.strip())
        print 'Read in '+str(len(mir_fcs))+' miRNA fold change values'

    #now load in miRNA targets (default to TS once I have new values)
    mirTargs={}
    try:
        mirTargs=pickle.load(open(opts.target_file,'rU'))
    except pickle.UnpicklingError:
        print '--miRNA-targets is not a proper pickle file, please try again'
        sys.exit()
    print 'Loaded up target data for '+str(len(mirTargs))+' families'

    #double check to make sure miRNA and targets have common names
    #print ','.join(miRs.keys())
    common=[a for a in miRs.keys() if a in mirTargs.keys()]
    if len(common)==0:
        print 'miRNA weights do not overlap with miRNA targets. This could be a naming strategy, please check files and try again'
        sys.exit()
    else:
        print 'Found '+str(len(common))+' miRNAs with target info'

    #now read in tf data
    tfs={}
    fn=opts.tf_file
    for row in open(fn,'rU').readlines():
        if len(row.strip().split('\t'))!=2:
            print '--TFWeights file needs to be tab-delimited with 2 columns, one for TF names, one for log fold change'
        
            sys.exit()
        else:
            m,val=row.strip().split('\t')
            tfs[m.strip()]=float(val.strip())
    print 'Found '+str(len(tfs))+' nodes that could be considered if they have binding sites in chromatin data...'

    #finally load up histones into dictionaries
    regions=opts.exp_names.split(',')
    up_regs=opts.upRegRegions_files.split(',')
    down_regs=opts.downRegRegions_files.split(',')
    if len(regions)!=len(up_regs) and len(regions)!=len(down_regs):
        print 'Error: --chromatinRegionNames needs to have same number of items as --upRegulatedRegions and --downRegulatedRegions'
        sys.exit()
    up_chroms,down_chroms={},{} ##dictionary of graphs!!
    for ind,val in enumerate(regions):
        print 'Loading '+val+' chromatin regions...'
        try:
            upg=pickle.load(open(up_regs[ind],'rU'))
            downg=pickle.load(open(down_regs[ind],'rU'))
        except pickle.UnpicklingError:
            'Print file found, but one of your '+val+' files cannot be upickled'
            sys.exit()

        up_chroms[val]=upg
        down_chroms[val]=downg
    
    upreg_genes=[a for a in mRNAs.keys() if mRNAs[a]>0]
    downreg_genes=[a for a in mRNAs.keys() if mRNAs[a]<0]

    ##first prepare mRNA inputs
    ##need to split up mRNAs by commodity, add mrna label
    conditions=['Up','Down']
    newWeights={}
    mirWeights={}
    absmirs={}
    upmirs,downmirs=[],[]
    if len(mir_fcs)==0: #then we treat original miRNA file as FC values
        upmirs=[a for a in miRs.keys() if miRs[a]>0]
        downmirs=[a for a in miRs.keys() if miRs[a]<0]
    else:
        upmirs=[a for a in mir_fcs.keys() if mir_fcs[a]>0]
        downmirs=[a for a in mir_fcs.keys() if mir_fcs[a]<0]
        
    for mir,val in miRs.iteritems():
        absmirs[mir]=np.abs(val)
    for c in conditions:
        newWeights[c]={}
        mirWeights[c]=absmirs
    for m in upreg_genes:
        newWeights['Up'][m+'mrna']=np.abs(mRNAs[m])
    for m in downreg_genes:
        newWeights['Down'][m+'mrna']=np.abs(mRNAs[m])
    
    if opts.do_ko:
        lo='TF'
    else:
        lo=''
        
    ##should I add this as command line option? or just use it by default since it will be more accurate?
    countflow=True
    
    print '----------------Building Network Graph-----------------------'
    graph,caps=build_samnet_network(miRs,mir_fcs,mirTargs,upreg_genes,downreg_genes,tfs,up_chroms,down_chroms,do_hier=True,cutoff=float(opts.histCutoff),draw=opts.plot_only)

    print '----------------Running SAMNet-----------------------'
    #now graph should be built, can actually run samnet!
    if not opts.plot_only:
        mirs,prots,newtfs,mrnas=runSamNetWithMirs(graph,newWeights,mirWeights,opts.gamma,opts.addpath,opts.out+'_gamma'+opts.gamma,conditions=conditions,leaveOut=lo,node_caps=caps,debug=False,sinkGamma=False)
        siffile=opts.out+'_gamma'+opts.gamma+'multiComm_mcfs.sif'
        #now make cytoscape output
        if opts.mirfc_file!='':
            mirclust=True
        else:
            mirclust=False
        make_cyto28_files(opts.mRNA_file,opts.miRNA_file,siffile,opts.out,mirclust)

    if opts.num_samps is not None:
        othermirs,otherprots,othertfs,othermrnas,totalreps=randomizeGraphAndRun(int(opts.num_samps),miRs=miRs, mirFC=mir_fcs,mirTargs=mirTargs,mirWeights=mirWeights,mrna_weights=newWeights,upreg_genes=upreg_genes,downreg_genes=downreg_genes,tfs=tfs,up_chroms=up_chroms,down_chroms=down_chroms,cutoff=float(opts.histCutoff),gamma=opts.gamma,addpath=opts.addpath,out=opts.out,conditions=conditions,orig_out=opts.out+'_gamma'+opts.gamma+'multiComm',levels=opts.levels.split(','),countFlow=countflow)
        #now we can compute stats!
        statsfile=open(opts.out+'gamma'+opts.gamma+'_'+opts.num_samps+'_randomization.xls','w')
        statsfile.write('Node Type\tNode\tFraction Selected\n')
        #first miRNAs
        for m in mirs.keys():
            if m in othermirs.keys():
                if countflow:
                    pval=float(len([a for a in othermirs[m] if a>mirs[m]]))/float(totalreps)
                else:
                    pval=float(othermirs[m])/float(totalreps)
                statsfile.write('miRNA\t%s\t%f\n'%(m,pval))
            else:
                statsfile.write('miRNA\t%s\t%f\n'%(m,float(0)))
        # for m in othermirs.keys():
        #     if m not in mirs:
        #         if countflow:
        #             statsfile.write('miRNA-randOnly\t%s\t%f\n'%(m,float(len(othermirs[m]))/float(totalreps)))
        #         else:
        #             statsfile.write('miRNA-randOnly\t%s\t%f\n'%(m,float(othermirs[m])/float(totalreps)))
        #then tfs
        for m in prots:
            if m in otherprots.keys():
                if countflow:
                    pval=float(len([a for a in otherprots[m] if a>prots[m]]))/float(totalreps)
                else:
                    pval=float(otherprots[m])/float(totalreps)
                statsfile.write('TF\t%s\t%f\n'%(m,pval))
            else:
                statsfile.write('TF\t%s\t%f\n'%(m,float(0)))
        # for m in otherprots.keys():
        #     if m not in prots:
        #         statsfile.write('TF-randOnly\t%s\t%f\n'%(m,float(otherprots[m])/float(totalreps)))
        #then motifs, I hope
        for m in newtfs:
            if m in othertfs.keys():
                if countflow:
                    pval=float(len([a for a in othertfs[m] if a>newtfs[m]]))/float(totalreps)
                else:
                    pval=float(othertfs[m])/float(totalreps)
                statsfile.write('Motif\t%s\t%f\n'%(m,pval))
            else:
                statsfile.write('Motif\t%s\t%f\n'%(m,float(0)))
        # for m in othertfs.keys():
        #     if m not in newtfs:
        #         statsfile.write('Motif-randOnly\t%s\t%f\n'%(m,float(othertfs[m])/float(totalreps)))
                
        statsfile.close()

if __name__=='__main__':
    main()
