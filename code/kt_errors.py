'''
python kt_errors.py dir/count dir/out metadata.csv experiment,replicate

'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import statsmodels.api as sm
import random


if len(sys.argv)!=6:
	raise Exception('Specify all files')

nboot=1000 #number of bootstrapped samples to get kappa std

countdir=sys.argv[1]
outdir=sys.argv[2]
meta=pd.read_csv(sys.argv[3])


exprep=sys.argv[4].split(',')
strain=exprep[0]
replicate=exprep[1][0]
freqs=pd.read_csv(countdir + '/{}_{}_freqs.csv'.format(strain,replicate))
Rtot=pd.read_csv(countdir+ '/{}_{}_Rtot.csv'.format(strain,replicate)).drop(columns=['Unnamed: 0']).set_index('Day').squeeze()

cfus=pd.read_csv(sys.argv[5])
print('Calculating noise model error parameters for {},{}'.format(strain,replicate))


cond1=meta['Strain']==strain
cond2=meta['Replicate']==np.int(replicate)
mm=meta[cond1 & cond2]


cond1=cfus['Strain']==strain
cond2=cfus['Replicate']==np.int(replicate)
cfus=cfus[cond1 & cond2]
cfus['avCFU']=cfus[['CFU1','CFU2']].mean(axis=1)/cfus['Dilution']
cfus2 = cfus[['avCFU','Day']].set_index('Day').squeeze()

#average cfus

# What are the adjacent days? 
# this is all very inefficient, sorrydaypairs=[]

days=np.sort(mm['Day'])
days = days[days>=0]
daypairs=[]
for i,d1 in enumerate(days[:-1]):
	d2=d1+1
	daypairs.append((d1,d2,np.int(d2)-np.int(d1)))

# Calculate kappa for each day pair
kdata={
	'Day 1':[],
	'Day 2':[],
	'kappa':[],
	'kappa std':[],
	'R1':[],
	'R2':[],
	'xbar':[],
	'xbar std':[],
	'ntransfers':[],
}

fmin={}
fmax={}
for day in days[:-1]:
	fmin[day]=np.max([30/Rtot[day], 30/cfus2[day] ])
	fmax[day]=7e-4 #500/Rtot[day]


for dp in daypairs:
	d1_i,d2_i,dt=dp
	d1=str(d1_i)
	d2=str(d2_i)
	R1=np.float(Rtot[d1_i])
	R2=np.float(Rtot[d2_i])
	cd=freqs[[d1,d2]]
	#print(cd)
	cond1=cd[d1]>fmin[d1_i]
	cond2=cd[d1]<fmax[d1_i]
	cd=cd[cond1 & cond2].sample(frac=1).reset_index(drop=True)

	kva=np.array(np.sqrt(cd[d1]) - np.sqrt(cd[d2]))

	xbar_v=(np.log(cd[d1]) - np.log(cd[d2]))/6.64
	xbar = np.median(xbar_v)

	kmed = np.median(kva)
	kappa_t = (np.median(np.abs(kva - kmed))/0.67449)**2
	kv=list(kva)
	boot_samples=[]
	xbar_samples=[]
	for j in range(nboot):
		kv_boot=np.array(random.choices(kv,k=len(kv)))
		kmed_boot = np.median(kv_boot)
		kappa_boot= (np.median(np.abs(kv_boot - kmed_boot))/0.67449)**2
		boot_samples.append(kappa_boot)

		x_boot=np.array(random.choices(xbar_v,k=len(xbar_v)))
		xbar_samples.append(np.median(x_boot))

	kdata['Day 1'].append(d1)
	kdata['Day 2'].append(d2)
	kdata['kappa'].append(kappa_t)
	kdata['kappa std'].append(np.std(boot_samples))
	kdata['xbar'].append(xbar)
	kdata['xbar std'].append(np.std(xbar_samples))
	kdata['R1'].append(R1)
	kdata['R2'].append(R2)
	kdata['ntransfers'].append(dt)


k_df=pd.DataFrame(kdata)
k_df.to_csv(outdir+'/{}_{}_kappa.csv'.format(strain,replicate))
print('Wrote '+outdir+'/{}_{}_kappa.csv'.format(strain,replicate))