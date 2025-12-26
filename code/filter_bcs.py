import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import statsmodels.api as sm
import random
import copy
import outliers_likelihoods as olh
import matplotlib.pyplot as plt

if len(sys.argv)!=6:
	raise Exception('Specify all files')


countdir=sys.argv[1]
meta=pd.read_csv(sys.argv[3])


exprep=sys.argv[4].split(',')
strain=exprep[0]
replicate=exprep[1][0]
counts=pd.read_csv(countdir + '/{}_{}_counts.csv'.format(strain,replicate)).set_index('super bc')
Rtot=pd.read_csv(countdir+ '/{}_{}_Rtot.csv'.format(strain,replicate)).drop(columns=['Unnamed: 0'],errors='ignore')
fitness=pd.read_csv(sys.argv[2] + '/{}_{}_fitness.csv'.format(strain,replicate))

try:
	counts.drop(columns=['-1'],inplace=True)
	Rtot=Rtot[Rtot['Day']>=0]
except:
	pass

qraw=0.5
qq_u=np.quantile(counts['0'],qraw)
qraw=0.1
qq_l=np.quantile(counts['0'],qraw)
counts=counts[counts['0']<qq_u]
counts=counts[counts['0']>qq_l]


fitness=fitness[fitness['super bc'].isin(list(counts.index))]

# filter out outliers
counts = counts.drop(index=list(fitness[fitness['pval']<5e-2]['super bc']))#.drop(columns='-1',errors='ignore')

cfus=pd.read_csv(sys.argv[5])
print('Filtering out outliers, high count bcs and merging low-count bcs for {},{}'.format(strain,replicate))


cond1=meta['Strain']==strain
cond2=meta['Replicate']==np.int(replicate)
mm=meta[cond1 & cond2]


cond1=cfus['Strain']==strain
cond2=cfus['Replicate']==np.int(replicate)
cfus=cfus[cond1 & cond2]
cfus['avCFU']=cfus[['CFU1','CFU2']].mean(axis=1)/cfus['Dilution']
cfus2 = cfus[['avCFU','Day']].set_index('Day').squeeze()

days=np.sort(mm['Day'])
days = days[days>=0]
days_s = [str(t) for t in days]



counts_arr = np.array(counts[days_s])

Rtot_new=np.sum(counts_arr,axis=0)

Rtot_old=np.array(Rtot['R'])

# get minimum counts such that CLT is likely to apply for both cell counts and sequence counts
fmin={}
mins=[]
for day in days[:-1]:
	fmin[str(day)]=np.max([50, Rtot_old[day]*50/cfus2[day] ])
	mins.append(Rtot_old[day]*50/cfus2[day])

fmin[str(days[-1])]=np.max([50,np.mean(mins)])

#print(freqs_arr[:,1])

for i,day in enumerate(days_s):
	while np.min(counts_arr[:,i])<fmin[day]:
		minloc = np.argmin(counts_arr[:,i])
		minarr = copy.deepcopy(counts_arr[minloc,:])
		counts_arr = np.delete(counts_arr,minloc,0)
		minloc2 = random.randint(0,len(counts_arr[:,i])-1)
		counts_arr[minloc2,:]+=minarr


freqs_arr=counts_arr/Rtot_old

d3=pd.DataFrame(freqs_arr,columns=days_s)
d3.to_csv(countdir+'/{}_{}_freqs_filtered_orig.csv'.format(strain,replicate),index=False)

freqs_arr=counts_arr/Rtot_new

d3=pd.DataFrame(freqs_arr,columns=days_s)
d3.to_csv(countdir+'/{}_{}_freqs_filtered.csv'.format(strain,replicate),index=False)

Rtot['R new']=Rtot_new
Rtot['R old/new']=Rtot['R']/Rtot['R new']

Rtot.to_csv(countdir+ '/{}_{}_Rtot.csv'.format(strain,replicate),index=False)


# add new Rtot to Rtot file here


#merge low-count barcodes into the next-lowest count bcs
'''
for day in days_s:
	while np.min(freqs[day])<fmin[day]:
		minloc=freqs[day].idxmin()
		minseries=copy.deepcopy(freqs.loc[minloc])
		freqs.drop(inplace=True,index=minloc)
		minloc2=freqs[day].idxmin()
		freqs.loc[minloc2]+=minseries


freqs.to_csv(countdir+'/{}_{}_freqs_filtered.csv'.format(strain,replicate))
'''