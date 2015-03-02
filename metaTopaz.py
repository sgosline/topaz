'''
metaTopaz.py module of the topaz mackage
Assembles  transcriptional data from garnet and gene expression data into one single network
to identify putative transcriptional regulations using correlation data

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
from collections import defaultdict

#import main topaz script
import topaz
global samnet

from itertools import tee, izip
def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def build_mirna_dis_network(miR,disTargs,disMrnaCors,disChromRegions,tfWeights={},do_hier=True,delim='.',cutoff=0.2):
    '''
    Constructs mRNA networks from miRNA_mRNA correlation values across various expression datasets
    to compare the impact of a single miRNA across multiple systems

    miR-TF edges are weighted by target scan prediction (if TF-miR have a PC<0)
    TF-histone edges are weighted by TF-miR anti-correlation
    histoine-mRNA are weighted by score
    mRNA-sink are weighted by anti-cor from miR in disease
    '''
    ##first start with all histone marks, and add in chrom-region-mRNA level of network
    hmarks=set()
    [hmarks.update(disChromRegions[d].keys()) for d in disChromRegions.keys()]
    
    print 'Evaluating '+str(len(hmarks))+' types of epigenetic data with '+str(len(disTargs))+' diseases'
    #distinguish TFs by whether or not they are up-regulated of down-regulated
    upperCaseTargs=defaultdict(dict)
    upperCaseCors=defaultdict(dict)
    for dis in disTargs.keys():
        geneval=disTargs[dis]
        for gene in geneval.keys():
            #corvals=disMrnaCors[dis]
            #if gene in corvals.keys():
            #    if corvals[gene]<0:
            upperCaseTargs[dis][gene.upper()]=geneval[gene]
            #    else:
            #        print miR+' and TF '+gene+' have cor value of '+corvals[gene]+' in '+dis+', skipping'
        geneval=disMrnaCors[dis]
        for gene in geneval.keys():
            upperCaseCors[dis][gene.upper()]=geneval[gene]
        print 'Have %d %s targets and %d anti-correlated genes for %s'%(len(upperCaseTargs[dis]),miR,len(upperCaseCors[dis]),dis)

    
    graph=nx.DiGraph()
    node_caps=defaultdict(dict)
    downTfs={}
    for d in upperCaseTargs.keys():
        #first collect list of anti-correlated TFs
        anticor=[tf for tf in upperCaseCors[d].keys() if upperCaseCors[d][tf]<0]
        #now filter those TF targets to be sure they are somewhat negatively correlated?
#        downTfs.update([tf for tf in upperCaseTargs[d].keys() if tf in anticor])
        newl=[tf for tf in upperCaseTargs[d].keys() if tf in anticor]
        downTfs[d]=newl

#        print 'Have '+str(len(newl))+' genes anti-correlated with '+miR+' in '+d
    #print 'Total of '+str(len(down_tfs))+' tfs:'+','.join(down_tfs)
    

    tfsWithBinding=set()
    #for each disease, we create a different transcriptional network, weighting
    #each tf by its anti-correlation
    for d in disTargs.keys():
        #iterate through dis commodity for every TF, select proper histone targets
        for hist in disChromRegions[d].keys():
            #iterate through every combination
            chromGraph=disChromRegions[d][hist]
            for matrix,targs in chromGraph.adjacency_iter():
                if len(targs)==0:
                    continue
                #first figure out which TFs are in the matrix
                tfs=[t.upper() for t in matrix.split(delim) if t in downTfs[d]]
                if len(tfs)==0:
                    continue
 #               print 'evaluating matrix: '+matrix+' with tfs '+','.join(tfs)
                #for each TF that is part of matrix, connect TF to chromatin node
                for tft in tfs:
                    if tft not in graph.nodes():
                        graph.add_edge('Dummy',tft)
                    if tft in upperCaseCors[d].keys():
                        tfw=np.abs(float(upperCaseCors[d][tft]))
                    graph.add_edge(tft,d+'.'+hist+'.'+matrix,weight=tfw)
                    tfsWithBinding.add(tft)
                for targ,eattr in targs.iteritems():
                    #for every target near chromatin region, add in anti-correlated miRNA
                    weight=eattr['weight']
                    if weight<0.5:
                        continue
#we don't need to check if there is anti-correlation data for target
                    #that will get filtered out in SAMNet
#                    mtarg=re.sub('mrna','',targ)
 #                   if mtarg in mirMrnaCorrelations[mir].keys() and mirMrnaCorrelations[mir][mtarg]<(-1*cutoff):
                    graph.add_edge(d+'.'+hist+'.'+matrix,d+'.'+targ,weight=weight)
        print 'Graph has '+str(graph.number_of_nodes())+' nodes and '+str(graph.number_of_edges())+' edges for '+str(len(tfsWithBinding))+' tfs with binding out of all down-regulated TFs'

    print 'Total selected tfs: '+','.join(tfsWithBinding)
    if do_hier:
        print 'Constructing node-specific capacties'
        ##Hierarchical capacities more fully mimic the flow model - miRNA-TF edges can handle more flow then TF-mRNA capacities...
        for node in graph.nodes():
            treat='Dummy'
            if node!=treat:
                spdist=1.0
                try:
                    spdist=nx.shortest_path_length(graph,treat,node,None)
                except nx.exception.NetworkXNoPath:
                    continue
                node_cap=np.power(10.0,float(spdist)*-1.0)
                #print treat,node,str(spdist),str(node_cap)
                for m in upperCaseCors.keys():
                    node_caps[m][node]=node_cap
    graph.remove_node('Dummy')

    #return network
    return graph,node_caps    

def build_mirna_cor_network(miRs,miRTargs,mirMrnaCorrelations,chromRegions,tfWeights={},do_hier=True,delim='.',cutoff=0.20):
    '''
    Constructs mRNA network from miRNA-mRNA correlation values in a single disease/systems
    miRNA-tf edges will be weighted by Target value (not in this function)
    TF-histone edges will be weighted by tfWeights dictionary - must not be commodity specific. If empty willd efault to 1
    histon-mRNA edges will be weighted by score
    mRNA-sink edges will be weighted by anti-correlation (act) or correlation (rep) with miRNA
    '''
    ##first start with all histone marks, and add in chrom-region-mRNA level of network
    hmarks=chromRegions.keys()
    
    print 'Evaluating '+str(len(hmarks))+' types of epigenetic data with '+str(len(miRs))+' miRNAs'
    #distinguish TFs by whether or not they are up-regulated of down-regulated
    upperCaseTargs=defaultdict(dict)
    for mir in miRs:
        geneval=miRTargs[mir]
        for gene in geneval.keys():
            upperCaseTargs[mir][gene.upper()]=geneval[gene]


    ##for each miRNA, get a list of mRNAs that are down-regulated below a specific cutoff
    down_tfs=set()


    graph=nx.DiGraph()
    node_caps=defaultdict(dict)
    for m in upperCaseTargs.keys():
        anticor=[tf.upper() for tf in mirMrnaCorrelations[m].keys() if mirMrnaCorrelations[m][tf]<(-1*cutoff)]
        newl=[tf for tf in upperCaseTargs[m].keys() if tf in anticor]
        down_tfs.update(newl)
        print 'Have '+str(len(newl))+' genes anti-correlated with '+m+' for threshold <-'+str(cutoff)
    #print 'Total of '+str(len(down_tfs))+' tfs:'+','.join(down_tfs)
    

    tfsWithBinding=set()
    #iterate through miRNA commodity for every TF, select proper histone targets
    for hist in hmarks:
        #iterate through every combination
        chromGraph=chromRegions[hist]
        for matrix,targs in chromGraph.adjacency_iter():
            if len(targs)==0:
                continue
            #first figure out which TFs are in the matrix
            tfs=[t.upper() for t in matrix.split(delim) if t in down_tfs]
            if len(tfs)==0:
                continue
            #print 'evaluating matrix: '+matrix+' with tfs '+','.join(tfs)
            #for each TF that is part of matrix, connect TF to chromatin node
            for tft in tfs:
                if tft not in graph.nodes():
                    graph.add_edge('Dummy',tft)
                if tft not in tfWeights.keys():
                    if len(tfWeights)==0:
                        tfw=1.0
                    else:
                        tfw=np.mean(tfWeights.values())
                else:
                    tfw=tfWeights[tft]
                graph.add_edge(tft,hist+'.'+matrix,weight=tfw)
                tfsWithBinding.add(tft)
            for targ,eattr in targs.iteritems():
                #for every target near chromatin region, add in anti-correlated miRNA
                weight=eattr['weight']
                if weight<0.5:
                    continue
                    #we don't need to check if there is anti-correlation data for target
                    #that will get filtered out in SAMNet
#                    mtarg=re.sub('mrna','',targ)
 #                   if mtarg in mirMrnaCorrelations[mir].keys() and mirMrnaCorrelations[mir][mtarg]<(-1*cutoff):
                graph.add_edge(hist+'.'+matrix,targ,weight=weight)
        print 'Graph has '+str(graph.number_of_nodes())+' nodes and '+str(graph.number_of_edges())+' edges for '+str(len(tfsWithBinding))+' tfs with binding out of '+str(len(down_tfs))+' down-regulated TFs'

    print 'Total selected tfs: '+','.join(tfsWithBinding)
    if do_hier:
        print 'Constructing node-specific capacties'
        ##Hierarchical capacities more fully mimic the flow model - miRNA-TF edges can handle more flow then TF-mRNA capacities...
        for node in graph.nodes():
            treat='Dummy'
            if node!=treat:
                spdist=1.0
                try:
                    spdist=nx.shortest_path_length(graph,treat,node,None)
                except nx.exception.NetworkXNoPath:
                    continue
                node_cap=np.power(10.0,float(spdist)*-1.0)
                #print treat,node,str(spdist),str(node_cap)
                for m in mirMrnaCorrelations.keys():
                    node_caps[m][node]=node_cap
    graph.remove_node('Dummy')

    #return network
    return graph,node_caps

#def run_meta_topaz(miRNAs,mirTargs,mirMrnaCor,chroms,gamma,path_to_samnet,outfile,tfWeights={},do_hier=False,do_ko=False,cutoff=0.25):
def run_meta_topaz(miRNAs,targScores,corScores,chroms,gamma,path_to_samnet,outfile,tfWeights={},do_hier=False,do_ko=False,cutoff=0.25,numiters=None):
    print '----------------Building Network Graph-----------------------'
    if len(miRNAs)==1:#we are building pan-disease network
        graph,caps=build_mirna_dis_network(miRNAs[0],targScores,corScores,chroms,do_hier=True,cutoff=cutoff)
    else:
        graph,caps=build_mirna_cor_network(miRNAs,targScores,corScores,chroms,do_hier=True,cutoff=cutoff)

    ##need to split up mRNAs by commodity, add mrna label
    # conditions=['Up','Down']
    # newWeights={}
    # mirWeights={}
    # absmirs={}
    # for mir,val in miRs.iteritems():
    #     absmirs[mir]=np.abs(val)
    # for c in conditions:
    #     newWeights[c]={}
    #     mirWeights[c]=absmirs
    # for m in upreg_genes:
    #     newWeights['Up'][m+'mrna']=np.abs(mRNAs[m])
    # for m in downreg_genes:
    #     newWeights['Down'][m+'mrna']=np.abs(mRNAs[m])
    newcors=defaultdict(dict)
    allmrnas=set()
    for mir,corval in corScores.iteritems():#this could be miRNA or disease specific
        for targ in corval.keys():
            if len(miRNAs)>1:
                newcors[mir][targ+'mrna']=corval[targ]
                allmrnas.add(targ+'mrna')
            else:
                newcors[mir][mir+'.'+targ+'mrna']=corval[targ]
                allmrnas.add(mir+'.'+targ+'mrna')
                
    print 'Pruning graph to remove spurious chromatin-mrna edges from '+str(graph.number_of_nodes())+' nodes and '+str(graph.number_of_edges())
    #now try to prune graph
    for n in graph.nodes():
        if 'mrna' in n and n not in allmrnas:
            graph.remove_node(n)
    print 'To '+str(graph.number_of_nodes())+' nodes and '+str(graph.number_of_edges())
            

    print '----------------Running SAMNet-----------------------'
    #now graph should be built, can actually run samnet!
    if do_ko:
        lo='TF'
        if numiters is not None:
            lo='NETWORK'
    else:
        if numiters is not None:
            lo='RNA'
        else:
            lo=''

    #this will create a standard run of SAMNet
    mirs,prots,tfs,mrnas=topaz.runSamNetWithMirs(network=graph,mir_weights=targScores,mrna_weights=newcors,gamma=gamma,outname=outfile+'_gamma'+gamma,conditions=targScores.keys(),samnet_path=path_to_samnet,leaveOut=lo,double=False,node_caps=caps,debug=False,sinkGamma=False,numiters=numiters)

    #if we do randomization, we need to sample the mrna_weights for each commodity, then re-run.
    #we do not want to save files for each randomization, we we create a temp directory and delete
    #we just want to keep a dictionary of how many times we select a TF?
    

def main():
    '''
    Parse arguments
    '''
    progdir=os.path.dirname(sys.argv[0])

    parser=OptionParser()
    ##collect mRNA expression data: tab-delimited file of differentially expressed mRNAs
    parser.add_option('--mRNAexp',dest='mRNA_file',type='string',help='Tab-delimited file containing 3 columns: miRNA commodity, anti-correlated mRNA and anti-correlation value')


    ##collect miRNA-Target data requires pickled dictionary of dictionaries
    parser.add_option('--miRNA-targets',dest='target_file',type='string',help='Tab-delimited or PKL file of dictionary of dictionaries containing miRNAs, targets and weights')
    
    ##which TFs should be considered?
    parser.add_option('--tfCorThreshold',dest='tfCor',type='string',help='Correlation threshold required to consider a TF as a regulator. Default is 0.25',default='0.25')

    ##needt o find universal tissue-specific weighting of TFS. For cancer I use median abundance - otherwise will be equal across TFs
    parser.add_option('--tfWeights',dest='tfweight',type='string',help='Tab-delimited file of TF names and weights. Otherwise all weights will default to 1.0',default='')
    
    ##collect TF-mRNA edges from GARNET output: data requires pickle files for each condition
    parser.add_option('--chromRegions',dest='chrom_regions_files',type='string',help='Comma-delimited list of TF-DNA networks (in PKL) for each histone or chromatin accessibility experiment for cell line')
    parser.add_option('--chromatinRegionNames',dest='exp_names',type='string',help='Comma-delimited list of names of chromatin experiments, e.g. H3K4me3')

###########################
    parser.add_option('--gamma',dest='gamma',type='string',help='SAMNet gamma parameter scales number of miRNAs selected')
    parser.add_option('--doTfKO',dest='do_ko',action='store_true',default=False,help='Do TF knockdown')
    parser.add_option('--outputPrefix',dest='out',type='string',help='Prefix for output',default='samnet')
    parser.add_option('--path-to-samnet',dest='addpath',type='string',default=os.path.join(progdir,'../SAMnet'),help='To run SAMNet we require path to SAMNet module')
    opts,args=parser.parse_args()

    ##first add samnet to path


    print '----------------Processing arguments-----------------------'
    #get mRNAs
    mf=open(opts.mRNA_file,'rU').readlines()
    mirMrnaCor=defaultdict(dict)
    for row in mf:
        if len(row.strip().split('\t'))!=3:
            print '--mRNAexp file needs to be tab-delimited with 3 columns, one for miRNA, one for mRNA, one forcorrelation value'
            sys.exit()
        else:
            mirna,mrna,val=row.strip().split('\t')
            mirMrnaCor[mirna][mrna]=float(val)
    print 'Read in '+str(len(mirMrnaCor))+' miRNAs from file: '+','.join(mirMrnaCor.keys())


    #now load in miRNA targets (default to TS once I have new values)
    mirTargs={}
    try:
        mirTargs=pickle.load(open(opts.target_file,'rU'))
    except pickle.UnpicklingError:
        print '--miRNA-targets is not a proper pickle file, trying to read as text'
        mirt=open(opts.target_file,'rU').readlines()
        for row in mirt:
            arr=row.strip().split('\t')
            if len(arr)!=3:
                print 'miR-target file needs to have 3 columns for miRNA, target and weight'
                sys.exit()
            mirTargs[arr[0].strip()]={arr[1].strip():float(arr[2].strip())}
            
    print 'Loaded up target data for '+str(len(mirTargs))+' families: '+','.join(mirTargs.keys())

    #double check to make sure miRNA and targets have common names
    common=[a for a in mirMrnaCor.keys() if a in mirTargs.keys()]
    if len(common)==0:
        print 'miRNA weights do not overlap with miRNA targets. This could be a naming strategy, please check files and try again'
        sys.exit()
    else:
        print 'Found '+str(len(common))+' miRNAs with target info'


    #load in tfweights
    tfWeights={}
    if opts.tfweight!='':
        for row in open(opts.tfweight,'rU').readlines():
            arr=row.strip().split('\t')
            if len(arr)!=2:
                print '--tfWeights needs to have 2 tab-delimited columns!'
                sys.exit()
            tfWeights[arr[0].strip()]=float(arr[1].strip())
        print 'Loaded in weights for '+str(len(tfWeights))+' TFS'
    
    #finally load up histones into dictionaries
    regions=opts.exp_names.split(',')
    regs=opts.chrom_regions_files.split(',')

    if len(regions)!=len(regs):
        print 'Error: --chromatinRegionNames needs to have same number of items as --chromRegions'
        sys.exit()
    chroms={}#dictionary of graphs
    for ind,val in enumerate(regions):
        print 'Loading '+val+' chromatin regions...'
        try:
            g=pickle.load(open(regs[ind],'rU'))
        except pickle.UnpicklingError:
            'Print file found, but one of your '+val+' files cannot be upickled'
            sys.exit()

        chroms[val]=g
        
    res=run_meta_topaz(miRNAs=common,mirTargs=mirTargs,mirMrnaCor=mirMrnaCor,chroms=chroms,gamma=opts.gamma,path_to_samnet=opts.addpath,outfile=opts.out,tfWeights=tfWeights,do_hier=True,do_ko=opts.do_ko,cutoff=opts.tfCor,numiters=None)
        
   # upreg_genes=[a for a in mRNAs.keys() if mRNAs[a]>0]
   # downreg_genes=[a for a in mRNAs.keys() if mRNAs[a]<0]


if __name__=='__main__':
    main()
