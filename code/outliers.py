
import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import time
import outliers_likelihoods as lh
from multiprocessing import Pool, TimeoutError
from functools import partial
from itertools import compress
import os
pd.set_option('mode.chained_assignment', None)

countdir=sys.argv[1]
kappadir=sys.argv[2]
outdir=sys.argv[3]
meta=pd.read_csv(sys.argv[4])
exprep=sys.argv[5].split(',')
strain=exprep[0]
replicate=exprep[1][0]
Rtot=pd.read_csv(countdir+ '/{}_{}_Rtot.csv'.format(strain,replicate)).drop(columns=['Unnamed: 0']).set_index('Day').squeeze()
counts=pd.read_csv(countdir+ '/{}_{}_counts.csv'.format(strain,replicate))
kappas=pd.read_csv(kappadir+ '/{}_{}_kappa.csv'.format(strain,replicate))


#select relevant metadata
cond1=meta['Strain']==strain
cond2=meta['Replicate']==np.int(replicate)
mm=meta[cond1 & cond2]

days=np.sort(mm['Day'])
tv = [str(t) for t in days[days>=0]]


#counts = counts[counts['gene_symbol']=='lacY']

print('Starting outlier detection for {},{}'.format(strain,replicate))

kappas['gens']=np.log(100)/np.log(2)



start=time.time()

def get_s_est(ggene):
	g,gene = ggene

	gene = gene[tv]

	if len(gene)<1:
		return None
	
	data_list=[]
	d00=lh.getgenedata(gene,kappas)
	if len(d00['r2'])>2:
		data_list.append(d00)

	if len(data_list)==0:
		return None
	s, s_std, pval_tot, pval_indiv, pval, success = lh.s_maxlikelihood(data_list)
	#print(s, s_std, s_pval, ll0)
	if success:
		return {
				'super bc':g,
				's':s,
				's std':s_std,
				'pval tot':pval_tot,
				'pval indiv':pval_indiv,
				'pval':pval,
				'stderr':s_std,
			}
	else:
		return None

batches = list(counts.groupby(by=['super bc']))
with Pool(processes=4) as pool:
	s_data = pool.map(get_s_est, batches)

print((time.time()-start)/3600)

s_data.append(None)
cond = pd.notnull(s_data)
s_data = list(compress(s_data, list(cond)))

df_save=pd.DataFrame(s_data)
s_signi,s_pval_corr=lh.FDR_correction(df_save['pval tot'])
df_save['corrected pval tot']=s_pval_corr

s_signi,s_pval_corr=lh.FDR_correction(df_save['pval indiv'])
df_save['corrected pval indiv']=s_pval_corr

s_signi,s_pval_corr=lh.FDR_correction(df_save['pval'])
df_save['corrected pval']=s_pval_corr

df_save.to_csv(outdir+'/{}_{}_fitness.csv'.format(strain,replicate))
print('Wrote '+outdir+'/{}_{}_fitness.csv'.format(strain,replicate))


#df_save[df_save['s significant']==1][['super bc']].to_csv(outdir+'/{}_{}_outliers.csv'.format(strain,replicate))

