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

import pickle,re,os,sys
import networkx as nx
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from copy import deepcopy

##need to update this if samnet is installed somewhere else
#sys.path.append("../SAMNet/src")



from itertools import tee, izip
def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def build_samnet_network(miRVals,miRTargs,upRegulatedGenes,downRegulatedGenes,tfsOfInterest,upRegulatedChromRegions,downRegulatedChromRegions,do_hier=True,delim='.',cutoff=0.5):
                         #top_mirs=15,upstream='',thresh=0.5,do_hier=True, default_node_cap=0.0,use_clip=True,use_cds_clip=False,dr_tfs={},ko_mrnas=[],wt_mrnas=[],ppi_file='',ppi_mapping_file='',hmarks=['H3K4me3','H3K27ac'],largeClusters=True,donorm=False):
    '''
    Assembles miRNAs and TFs and histone data into network to enable running by SAMNet
    miRTargs: dictionary of dictionary of miRNAs, targets and scores
    upRegulatedGenes: list of up-regulated genes
    downRegulatedGenes: list of down-regulated genes.
    tfsOfInterest: list of transcription factors of interest
    upRegulatedChromRegions: dictionary of TF-mRNA networks. Keys indiciate histone mark
    downRegulatedChromRegions: dictionary of TF-mRNA networks. Keys indiciate histone mark
    do_hier: add in hierarchical edge capacities
    delim: delimiter used in GARNET To collapse TF names
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
    up_mirs=[a for a in miRVals.keys() if miRVals[a]>0]
    down_mirs=[a for a in miRVals.keys() if miRVals[a]<0]
    print 'Have '+str(len(up_mirs))+' up-regulated miRs and '+str(len(down_mirs))+' down-regulated'
    #assemble weighted graph
    graph=nx.DiGraph()

    #now start with miRNas
    tfs_in_play=set()
    for m in miRTargs.keys():
        targs=miRTargs[m].keys()
        ##allow up-regulation and down-regulation to be explained by miRNAs in either direction
        graph.add_edge('Up',m,weight=1.0)###remove this node at end, just a dummy for sp computation
        graph.add_edge('Down',m,weight=1.0)###remove this node at end, just a dummy for sp computation
        for targ in targs:
            if m in up_mirs and targ.upper() in down_tfs:
                graph.add_edge(m,targ.upper(),weight=miRTargs[m][targ])#/max_vals)
                tfs_in_play.add(targ.upper())
            elif m in down_mirs and targ.upper() in up_tfs:
                graph.add_edge(m,targ.upper(),weight=miRTargs[m][targ])#/max_vals)
                tfs_in_play.add(targ.upper())
#            else:
#                print 'No edge between '+m+' and '+targ+' because they are not expressed in opposite directions'

    print 'Graph has '+str(len(graph.nodes()))+' nodes and '+str(len(graph.edges()))+' edges from '+str(len(miRTargs))+' miRNAs to '+str(len(tfs_in_play))+' possible TFs'

    ##then create TF-cond.Tf.histone interactions
    ###first start with Up######################################
    
    cond='Up'
    tfs_with_binding=set()
    for hist in upRegulatedChromRegions.keys():
        #for every histone create individual node for each matrix
        chromGraph=upRegulatedChromRegions[hist]
        for matrix,targs in chromGraph.adjacency_iter():
            #figure out which TFs target this matrix
            tfs=[t for t in matrix.split(delim) if t in tfs_in_play]
            ##add edge from every TF to the matrix
            for a in tfs:
                weight=np.abs(upperCaseTfs[a])##take absolute value of fold change, relying on miRNAs to filter
#                print 'Weight from '+a+' to '+matrix+' is :'+str(weight)
                if weight>0.0:
                    graph.add_edge(a,hist+'.'+cond+'.'+matrix,weight=weight)
                tfs_with_binding.add(a)
            #now iterate through every target and add edge from matrix to target
            for targ,eattr in targs.items():
                weight=eattr['weight']
                if weight>cutoff and targ in upRegulatedGenes:
                    graph.add_edge(hist+'.'+cond+'.'+matrix,targ+'mrna',{'weight':weight})                    
        print 'Graph has '+str(len(graph.nodes()))+' and '+str(len(graph.edges()))+' after '+hist+' Up regions and mrnas added'

    cond='Down'
    for hist in downRegulatedChromRegions.keys():
        #for every histone create individual node for each matrix
        chromGraph=downRegulatedChromRegions[hist]
        for matrix,targs in chromGraph.adjacency_iter():
            #figure out which TFs target this matrix
            tfs=[t for t in matrix.split(delim) if t in tfs_in_play]
            ##add edge from every TF to the matrix
            for a in tfs:
                weight=np.abs(upperCaseTfs[a])##take absolute value of fold change, relying on miRNAs to filter
                if weight>0.0:
                    graph.add_edge(a,hist+'.'+cond+'.'+matrix,weight=weight)
                tfs_with_binding.add(a)
            #now iterate through every target and add edge from matrix to target
            for targ,eattr in targs.items():
                weight=eattr['weight']
                if weight>cutoff and targ in downRegulatedGenes:
                    graph.add_edge(hist+'.'+cond+'.'+matrix,targ+'mrna',{'weight':weight})                    
        print 'Graph has '+str(len(graph.nodes()))+' and '+str(len(graph.edges()))+' after '+hist+' Down regions and mrnas added'
    print 'Only '+str(len(tfs_with_binding))+' tfs have binding out of '+str(len(tfs_in_play))+' mRNA above activity threshold'
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

    if 'Up' in graph.nodes():
        graph.remove_node('Up')
    if 'Down' in graph.nodes():
        graph.remove_node('Down')

    return graph,node_caps



def runSamNetWithMirs(network,mrna_weights,mir_weights,gamma,outname,conditions=['Up','Down'],leaveOut='',double=False,node_caps={},debug=False,sinkGamma=False):
    '''
    runs samnet with mirs in network - by this point mRNA nodes should have 'mrna' appended to their name to prevent direct edges...
    '''

    ##this calls the original SAMNet code - make sure you know where it is!
    flow,phens,prots,tfs,mrnas=samnet.run_rn(PPI_with_weights=network,indirect_weights=mir_weights,direct_weights={},graph_tr=nx.DiGraph(),mrna_weights=mrna_weights,output=outname,updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=debug,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)

    ##phens is another word for miRNAs - keep track of everythign
    total_phens=phens
    total_tfs=tfs
    total_mrnas=mrnas
    if len(prots)==0:
        res_tfs=[t for t in tfs.difference(phens)]
    else:        
        res_tfs=[t for t in prots.difference(tfs)]
  #  print 'Orig tfs '+','.join(tfs)
    print 'Res tfs '+','.join(res_tfs)


    ##collect list of output file prefixes to combine...
    outlist=[outname+'multiComm']
    mul_c='multiComm'
    ##This code will leave out single or double TFs at once
    if leaveOut.lower()=='tf':
#        rtfs=[t for t in tfs.difference(phens)]
        if double: #leave out 2 tfs at once
            for tfpair in pairwise(res_tfs):
                print 'running with '+tfpair[0]+' and '+tfpair[1]+' removed'
#                if not os.path.exists(outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED'+mul_c+'_ppi_attributes.eda'):
                newnetwork=deepcopy(network)
                newnetwork.remove_node(tfpair[0])
                newnetwork.remove_node(tfpair[1])                    
                newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr={},mrna_weights=mrna_weights,output=outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
                total_phens=total_phens.union(newphens)
                total_tfs=total_tfs.union(newtfs)
                total_mrnas=total_mrnas.union(newmrnas)
                outlist.append(outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED'+mul_c+'_edge_type.eda')
                
        else: #leave out 1 TF at once - this is what I do in the paper            
            for tf in res_tfs:
                print 'running with '+tf+' removed'
                #if not os.path.exists(outname+'_'+tf+'_REMOVED'+mul_c+'_ppi_attributes.eda'):
                newnetwork=deepcopy(network)
                newnetwork.remove_node(tf)
                newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr={},mrna_weights=mrna_weights,output=outname+'_'+tf+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
                total_phens=total_phens.union(newphens)
                total_tfs=total_tfs.union(newtfs)
                total_mrnas=total_mrnas.union(newmrnas)
                outlist.append(outname+'_'+tf+'_REMOVED'+mul_c+'_edge_type.eda')

    elif leaveOut=='mir':#code can also remove miRNAs from network, haven't really explored this
        if double:
            for tfpair in pairwise(phens):
                print 'running with '+tfpair[0]+' and '+tfpair[1]+' removed'
                newnetwork=deepcopy(network)
                newnetwork.remove_node(tfpair[0])
                newnetwork.remove_node(tfpair[1])
                newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr={},mrna_weights=mrna_weights,output=outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
#                outlist.append(outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED'+mul_c+'_ppi_attributes.eda')
                outlist.append(outname+'_'+tfpair[0]+'and'+tfpair[1]+'_REMOVED'+mul_c+'_edge_type.eda')
                total_phens=total_phens.union(newphens)
                total_tfs=total_tfs.union(newtfs)
                total_mrnas=total_mrnas.union(newmrnas)
                print newphens
        else:
            for mir in phens:
                newnetwork=deepcopy(network)
                newnetwork.remove_node(mir)
                newflow,newphens,newprots,newtfs,newmrnas=samnet.run_rn(PPI_with_weights=newnetwork,indirect_weights=mir_weights,direct_weights={},graph_tr={},mrna_weights=mrna_weights,output=outname+'_'+mir+'_REMOVED',updateIds='none',cap=1.0,gamma=gamma,solver='cplexamp',noTfs=True,debug=False,treat_weights=dict(),diff_ex_vals=dict(),de_cap='sink',doMCF=True,makeCombined=False,node_caps=node_caps,sinkGamma=sinkGamma)
                #                outlist.append(outname+'_'+mir+'_REMOVED'+mul_c+'_ppi_attributes.eda')
                total_phens=total_phens.union(newphens)
                total_tfs=total_tfs.union(newtfs)
                total_mrnas=total_mrnas.union(newmrnas)
                outlist.append(outname+'_'+mir+'_REMOVED'+mul_c+'_edge_type.eda')
                print newphens
            
    #now combine outlist files into a single unified sif file
    if len(outlist)>1:
        combined=samnet.combine_single_flows_to_make_multi(outlist[1:],outlist[0],collapse_edges=True,ismcf=True)
    return total_phens,prots,total_tfs,total_mrnas

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
    ##TODO: enable negative values for miRNAs to correspond to TFs that are increased in value
    
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
    parser.add_option('--doTfKO',dest='do_ko',action='store_true',default=False,help='Do TF knockdown')
    parser.add_option('--outputPrefix',dest='out',type='string',help='Prefix for output',default='samnet')
    parser.add_option('--path-to-samnet',dest='addpath',type='string',default=os.path.join(progdir,'../SAMnet'),help='To run SAMNet we require path to SAMNet module')
    opts,args=parser.parse_args()

    ##first add samnet to path
    sys.path.insert(0,opts.addpath)
    sys.path.insert(0,opts.addpath+'src')
    import samnet  ##call run_samnet command manually

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
    mif=open(opts.miRNA_file,'rU')
    for row in mif:
        if len(row.strip().split('\t'))!=2:
            print '--miRNAs file needs to be tab-delimited with 2 columns, one for miR one for lfc'
            sys.exit()
        else:
            mir,val=row.strip().split('\t')
            miRs[mir]=float(val)
    print 'Read in '+str(len(miRs))+' miRNA values'

    #now load in miRNA targets (default to TS once I have new values)
    mirTargs={}
    try:
        mirTargs=pickle.load(open(opts.target_file,'rU'))
    except pickle.UnpicklingError:
        print '--miRNA-targets is not a proper pickle file, please try again'
        sys.exit()
    print 'Loaded up target data for '+str(len(mirTargs))+' families'

    #double check to make sure miRNA and targets have common names
    common=[a for a in miRs.keys() if a in mirTargs.keys()]
    if len(common)==0:
        print 'miRNA weights do not overlap with miRNA targets. This could be a naming strategy, please check files and try again'
        sys.exit()
    else:
        print 'Found '+str(len(common))+' miRNAs with target info'

    #now read in tf data
    tfs={}
    for row in open(opts.tf_file,'rU').readlines():
        if len(row.strip().split('\t'))!=2:
            print '--TFWeights file needs to be tab-delimited with 2 columns, one for TF names, one for log fold change'
        
            sys.exit()
        else:
            m,val=row.strip().split('\t')
            tfs[m]=float(val)
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

    print '----------------Building Network Graph-----------------------'
    graph,caps=build_samnet_network(miRs,mirTargs,upreg_genes,downreg_genes,tfs,up_chroms,down_chroms,do_hier=True)

    ##need to split up mRNAs by commodity, add mrna label
    conditions=['Up','Down']
    newWeights={}
    mirWeights={}
    for c in conditions:
        newWeights[c]={}
        mirWeights[c]=miRs
    for m in upreg_genes:
        newWeights['Up'][m+'mrna']=np.abs(mRNAs[m])
    for m in downreg_genes:
        newWeights['Down'][m+'mrna']=np.abs(mRNAs[m])
    

    print '----------------Running SAMNet-----------------------'
    #now graph should be built, can actually run samnet!
    if opts.do_ko:
        lo='TF'
    else:
        lo=''
    mirs,prots,tfs,mrnas=runSamNetWithMirs(graph,newWeights,mirWeights,opts.gamma,opts.out+'_gamma'+opts.gamma,conditions=conditions,leaveOut=lo,double=False,node_caps=caps,debug=False,sinkGamma=False)

if __name__=='__main__':
    main()
