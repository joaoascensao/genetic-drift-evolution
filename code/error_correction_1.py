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
from polyleven import levenshtein
import Levenshtein

if len(sys.argv)!=5:
	raise Exception('Specify codes and close files')

def getGC(seq):
    return (sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / 20)

gc_cutoff=0.2

codes_dir=sys.argv[1]
outdir=sys.argv[2]
meta=pd.read_csv(sys.argv[3])

exprep=sys.argv[4].split(',')
strain=exprep[0]
replicate=np.int(exprep[1][0])

meta=meta[(meta['Strain']==strain) & (meta['Replicate']==replicate)]

print('Starting error correction 1 for {},{}'.format(strain, replicate))



def write_list(ll,primer,filename):
	with open(outdir+'/{}_'.format(primer)+filename+'.txt', 'w') as f:
		for item in ll:
			f.write("{}\n".format(item))

def write_list2(ll,primer,filename):
	with open(outdir+'/{}_'.format(primer)+filename+'.txt', 'w') as f:
		for i1,i2 in ll:
			f.write("{},{}\n".format(i1,i2))

def write_list3(ll,strain,replicate,filename):
	with open(outdir+'/{}_{}_'.format(strain,replicate)+filename+'.txt', 'w') as f:
		for item in ll:
			f.write("{}\n".format(item))

all_bc_nodes=[]

for i,row in meta.iterrows():
	bc_nodes=[]
	bc_edges=[]

	high_count_nodes=[]

	num_primer=int(row['Primer'])
	if int(num_primer)<10:
		primer='IT00'+str(int(num_primer))
	else:
		primer='IT0'+str(int(num_primer))
	#primer_name='IT0'+str(int(num_primer))
	
	codes_df = pd.read_csv(codes_dir+'/{}.codes'.format(primer),sep='\t')
	close = pd.read_csv(codes_dir+'/{}.close'.format(primer),sep='\t')
	close1=list(close['code1'])
	close2=list(close['code2'])

	bc_nodes+=list(codes_df['barcode'])
	high_count_nodes+=list(codes_df[codes_df[primer]>=20]['barcode'])
	bc_edges+=list(zip(close1,close2))
	all_bc_nodes+=list(codes_df['barcode'])
		

	bc_nodes = list(set(bc_nodes))
	high_count_nodes = list(set(high_count_nodes))


	# exclude barcodes with high and low gc content
	gc_all = {seq:getGC(seq) for seq in bc_nodes}
	bc_nodes = [seq for seq in bc_nodes if (gc_all[seq]<1-gc_cutoff and gc_all[seq]>gc_cutoff)]
	high_count_nodes = [seq for seq in high_count_nodes if (gc_all[seq]<1-gc_cutoff and gc_all[seq]>gc_cutoff)]
	bc_edges = [(seq1,seq2) for seq1,seq2 in bc_edges if (gc_all[seq1]<1-gc_cutoff and gc_all[seq1]>gc_cutoff and gc_all[seq2]<1-gc_cutoff and gc_all[seq2]>gc_cutoff)]

	G = nx.Graph()
	G.add_nodes_from(bc_nodes)
	G.add_edges_from(bc_edges)


	#print(len([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]))
	#print(len(bc_nodes),len(bc_edges),len(high_count_nodes))


	ccs=[list(c) for c in nx.connected_components(G)]

	cc0_l=[]

	gg=0
	for cc in nx.connected_components(G):
		if len(cc)==1:
			cc0_l.append(list(cc)[0])
			gg+=1
	print(num_primer,gg,len(high_count_nodes))

	write_list(cc0_l,num_primer,'cc0')
	write_list(high_count_nodes,num_primer,'high_count_nodes')
	write_list2(bc_edges,num_primer,'edges0')
	del G

gc_all = {seq:getGC(seq) for seq in all_bc_nodes}
all_bc_nodes = list(set(all_bc_nodes))
all_bc_nodes = [seq for seq in all_bc_nodes if (gc_all[seq]<1-gc_cutoff and gc_all[seq]>gc_cutoff)]
write_list3(bc_nodes,strain,replicate,'all_nodes')

print('Done with error correction 1 for {},{}'.format(strain, replicate))


