'''

Usage: python codes_dir outdir meta exprep
'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys

if len(sys.argv)!=6:
	raise Exception('Specify codes and close files')



codes_dir=sys.argv[1]
errcorr_dir=sys.argv[2]
outdir=sys.argv[3]
meta=pd.read_csv(sys.argv[4])

exprep=sys.argv[5].split(',')
strain=exprep[0]
replicate=np.int(exprep[1][0])

meta=meta[(meta['Strain']==strain) & (meta['Replicate']==replicate)]

superbcs = pd.read_csv(errcorr_dir+'/{}_{}_superbarcodes.csv'.format(strain, replicate))

print('Counting barcodes for {},{}'.format(strain, replicate))

superbcs_codes=list(superbcs['super bc'])

days=[]
R=[]


def find_code(cd,x):
	if x in cd:
		return cd[x]
	else:
		return 0


for i,row in meta.iterrows():
	num_primer=int(row['Primer'])
	day=int(row['Day'])

	if int(num_primer)<10:
		primer='IT00'+str(int(num_primer))
	else:
		primer='IT0'+str(int(num_primer))
	#primer_name='IT0'+str(int(num_primer))
	codes_df = pd.read_csv(codes_dir+'/{}.codes'.format(primer),sep='\t')

	days.append(day)
	R.append(np.sum(codes_df[primer]))

	codes_dic=codes_df.set_index('barcode').squeeze().to_dict()
	#print(codes_dic)
	superbcs[day]=superbcs['barcode'].apply(lambda x: find_code(codes_dic,x))

superbcs.drop(columns=['barcode','av GC','Unnamed: 0'],inplace=True)
r=superbcs.groupby(by='super bc').sum().reset_index()


r['drop']=np.sum(np.sign( np.clip(r[days] - 15,0,1) ),axis=1)

r2=r[r['drop']>3].drop(columns=['drop']).reset_index().drop(columns=['index'])
r2.to_csv(outdir+'/{}_{}_counts.csv'.format(strain, replicate))

f=r2.copy()
f[days]=r2[days]/R
f.to_csv(outdir+'/{}_{}_freqs.csv'.format(strain, replicate))

#pd.Series(R,index=days)
pd.DataFrame({'Day':days,'R':R}).to_csv(outdir+'/{}_{}_Rtot.csv'.format(strain, replicate))

print('Done counting barcodes for {},{}'.format(strain, replicate))



