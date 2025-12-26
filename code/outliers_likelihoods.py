
import pandas as pd
import scipy as sp 
import numpy as np 
from scipy import optimize
from scipy.optimize import minimize
from scipy import interpolate
import scipy.special
import copy
import re
import random
from statsmodels.stats import multitest


largenum = 1e9

# clean up all data for a gene
def getgenedata(countdf,kappas):
	gene_data=kappas.copy()
	r1=[]
	r2=[]
	for i,row in kappas.iterrows():
		r1.append(np.float(countdf[str(int(row['Day 1']))]))
		r2.append(np.float(countdf[str(int(row['Day 2']))]))

	gene_data['r1']=r1
	gene_data['r2']=r2
	gene_data['kappa_t'] = 4*gene_data['kappa']*gene_data['R2']
	gene_data['r1_t'] = gene_data['r1']*gene_data['R2']/gene_data['R1']

	gene_data_dic = gene_data[['r2','kappa_t','r1_t','xbar']].to_dict(orient='series')
	
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))

	return gene_data_dic2



# stirlings approx
def gamma_stirlings(x):
	return 0.5*np.log(2*np.pi) + (x - 0.5)*np.log(x) - x + 1/(12*x)

# likelihood for a single barcode, single timepoint
def nb(y,mu,k):
	mu=mu+1e-2
	phi = mu/(k-1)
	g1 = y+phi
	g2 = y+1
	g3 = phi
	if np.all(g1>1): #stirlings approx
		dg1 = gamma_stirlings(g1)
	else:
		dg1 = np.longdouble(sp.special.loggamma(np.float64(g1)))
	
	if g2>1: #stirlings approx
		dg2 = gamma_stirlings(g2)
	else:
		dg2 = np.longdouble(sp.special.loggamma(np.float64(g2)))

	if np.all(g3)>1: #stirlings approx
		dg3 = gamma_stirlings(g3)
	else:
		dg3 = np.longdouble(sp.special.loggamma(np.float64(g3)))
	

	return dg1 - dg2 - dg3 + y*np.log(mu) - np.log(mu+phi)*(y+phi) + phi*np.log(phi)


def norm_pdf(y,mu,k):
	v = k*mu
	return -0.5*((y-mu)**2)/v - 0.5*np.log(2*np.pi*v)

# add likelihoods over all barcodes within a gene
def add_likelihoods(s,data_list):
	data = data_list[0]
	s=np.longdouble(s)

	#print(shareddata,data_list)
	
	#print(data)
	ll=np.zeros_like(s)
	ll_indv=[]
	for i in range(len(data['kappa_t'])):
		nb00=nb(data['r2'][i],data['r1_t'][i]*np.exp((s-data['xbar'][i])*6.64),data['kappa_t'][i])
		ll+=nb00
		ll_indv.append(nb00)
		
	return ll, ll_indv


# p-value calculation
def posterior_pval(svec,ll):
	
	log_post_int= interpolate.interp1d(svec,ll,kind='cubic')
	svec1 = np.linspace(min(svec),max(svec),5*10**4)
	ll00 = log_post_int(0)
	log_post2 = log_post_int(svec1)

	lam = np.nan_to_num(ll00-log_post2)
	mm=np.max(log_post2)
	post2=np.nan_to_num(np.exp(log_post2-mm))
	post3=post2/np.sum(post2)

	pval=np.sum(post3[lam>0])


	amx=np.argmax(log_post2)
	llm=np.max(log_post2)
	s_hat = svec1[amx]
	if np.abs(s_hat)==0.3:
		return np.nan, np.nan,np.nan, np.nan, False
	else:
		success=True

	# get 68% CI
	pvals=[]
	for si in svec:
		ll0 = log_post_int(si)
		lam = np.nan_to_num(ll0-log_post2)
		mm=np.max(log_post2)
		post2=np.nan_to_num(np.exp(log_post2-mm))
		post3=post2/np.sum(post2)
		pvals.append(np.sum(post3[lam>0]) - 0.3174)
	

	w=np.where(np.diff(np.sign(pvals))!=0)[0]
	if len(w)!=2:
		sl=np.nan
		su=np.nan
	else:
		sl=svec[w[0]]
		su=svec[w[1]]

	return s_hat, pval, sl, su, success

def posterior_pval2(svec,ll):
	
	log_post_int= interpolate.interp1d(svec,ll,kind='cubic')
	svec1 = np.linspace(min(svec),max(svec),5*10**4)
	ll00 = log_post_int(0)
	log_post2 = log_post_int(svec1)

	lam = np.nan_to_num(ll00-log_post2)
	mm=np.max(log_post2)
	post2=np.nan_to_num(np.exp(log_post2-mm))
	post3=post2/np.sum(post2)

	pval=np.sum(post3[lam>0])

	return pval

# get maximum likelihood estimate of fitness
def s_maxlikelihood(data_list, boot=False):

	s_v=np.linspace(-0.3,0.3,10**3)
	ll,ll_indv = add_likelihoods(s_v,data_list)
	
	s_hat, pval_tot,sl,su,success=posterior_pval(s_v,ll)

	pvals=[]
	for lli in ll_indv:
		pvals.append(posterior_pval2(s_v,lli))
	#print(pvals)

	return s_hat, su-sl, pval_tot, np.min(pvals), np.min(pvals + [pval_tot]), success
	

# BH FDR correaction
def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr
