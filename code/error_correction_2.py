'''

Usage: python codes_dir outdir meta exprep
'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import networkx as nx
import time
from multiprocessing import Pool, TimeoutError
from functools import partial
from itertools import compress


if len(sys.argv)!=5:
	raise Exception('Specify codes and close files')


cc_cutoff=1000

codes_dir=sys.argv[1]
outdir=sys.argv[2]
meta=pd.read_csv(sys.argv[3])

exprep=sys.argv[4].split(',')
strain=exprep[0]
replicate=np.int(exprep[1][0])

meta=meta[(meta['Strain']==strain) & (meta['Replicate']==replicate)]


print('Starting error correction 2 for {},{}'.format(strain, replicate))


allnodes = list(pd.read_csv(outdir+'/{}_{}_all_nodes'.format(strain,replicate)+'.txt',names=['bc'], header=None)['bc'])
alledges=[]
for i,row in meta.iterrows():
	num_primer=int(row['Primer'])

	edges1 = pd.read_csv(outdir+'/{}_edges0'.format(num_primer)+'.txt',names=['bc1','bc2'], header=None)
	edges2 = pd.read_csv(outdir+'/{}_leven_out'.format(num_primer)+'.txt',names=['bc1','bc2','dist'], header=None)
	edges1['dist']=1
	alledges.append(edges1)
	alledges.append(edges2)

edges=pd.concat(alledges)
edges=edges[edges['dist']<=3] #only include edges with lev dist<=3
e1=list(edges['bc1'])
e2=list(edges['bc2'])

G = nx.Graph()
G.add_nodes_from(allnodes)
G.add_edges_from(list(zip(e1,e2)))

ccs=[c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
del G

'''
G3 = nx.Graph()
G3.add_nodes_from(allnodes)

e1=list(edges['bc1'])
e2=list(edges['bc2'])
G3.add_edges_from(list(zip(e1,e2)))





print('largest cc:',np.max([len(c) for c in ccs]))

ccs1=[]
for cc in ccs:
	if len(cc)>cc_cutoff:
		S = G3.subgraph(cc).copy()
		cc_temp=[c for c in sorted(nx.connected_components(S), key=len, reverse=True)]
		ccs1+=cc_temp
		del S
	else:
		ccs1.append(cc)
del G3

print('largest cc:',np.max([len(c) for c in ccs1]))

#only include edges with lev dist<=2
G2 = nx.Graph()
G2.add_nodes_from(allnodes)
edges=edges[edges['dist']<=2]
e1=list(edges['bc1'])
e2=list(edges['bc2'])
G2.add_edges_from(list(zip(e1,e2)))

ccs2=[]
for cc in ccs1:
	if len(cc)>cc_cutoff:
		S = G2.subgraph(cc).copy()
		cc_temp=[c for c in sorted(nx.connected_components(S), key=len, reverse=True)]
		ccs2+=cc_temp
		del S
	else:
		ccs2.append(cc)

del G2

print('largest cc:',np.max([len(c) for c in ccs2]))


#only include edges with lev dist<=1
G1 = nx.Graph()
G1.add_nodes_from(allnodes)
edges=edges[edges['dist']<=1]
e1=list(edges['bc1'])
e2=list(edges['bc2'])
G1.add_edges_from(list(zip(e1,e2)))


ccs3=[]
for cc in ccs2:
	if len(cc)>cc_cutoff:
		S = G1.subgraph(cc).copy()
		cc_temp=[c for c in sorted(nx.connected_components(S), key=len, reverse=True)]
		ccs3+=cc_temp
		del S
	else:
		ccs3.append(cc)

del G1

print('largest cc:',np.max([len(c) for c in ccs3]))

#https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.stoerwagner.stoer_wagner.html#networkx.algorithms.connectivity.stoerwagner.stoer_wagner
#https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html#networkx.algorithms.community.modularity_max.greedy_modularity_communities
'''


def getGC(seqs):
    return [(sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / 20) for seq in seqs]

i=0
bcs=[]
sbc=[]
gc_cont=[]
for cc in ccs:
	bcs+=cc
	sbc+=[i]*len(cc)
	gc_cont+=getGC(cc)
	i+=1


pd.DataFrame({'super bc':sbc, 'barcode':bcs, 'av GC':gc_cont}).to_csv(outdir+'/{}_{}_superbarcodes.csv'.format(strain,replicate))

print('Done with error correction 2 for {},{}'.format(strain, replicate))
