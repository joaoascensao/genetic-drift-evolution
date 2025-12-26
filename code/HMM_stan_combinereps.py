import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pystan
import sys


stancode = """
functions {
  real LPTN_lpdf(real x, real mu, real std, real xi) {
    return normal_lpdf( x | mu, std)*xi;
  }
}

data {
  int<lower=0> T;           // # time points (equally spaced)
  int<lower=0> M_rep1;           // # of barcodes
  real<lower=0, upper=1> xi_rep1;
  vector[T] f_rep1[M_rep1];     // sqrt-transformed observed frequencies
  vector[T] R_rep1;              // total sequence reads for each time point
  vector[T-1] D_rep1;              // cfu dilution rate for each time point, after accounting for the transfer dilution
  int cfus1_rep1[T-1];           // CFUs for each time point (not multiplied by dilution rate)
  int cfus2_rep1[T-1];           // CFUs for each time point (not multiplied by dilution rate)

  int<lower=0> M_rep2;           // # of barcodes
  real<lower=0, upper=1> xi_rep2;
  vector[T] f_rep2[M_rep2]; // sqrt-transformed observed frequencies
  vector[T] R_rep2;              // total sequence reads for each time point
  vector[T-1] D_rep2;              // cfu dilution rate for each time point, after accounting for the transfer dilution
  int cfus1_rep2[T-1];           // CFUs for each time point (not multiplied by dilution rate)
  int cfus2_rep2[T-1];           // CFUs for each time point (not multiplied by dilution rate)
}
parameters {
  real<lower=1e-4,upper=50> sig;               // sigma^2, descendant number variance
  real<lower=1e-4,upper=30> eta_cfus; // error parameter for CFUs - 1

  vector<lower=1,upper=50>[T] c_rep1;               // c at each time point
  vector<lower=1e5,upper=5e8>[T-1] Nb_rep1;       // true bottleneck size at each time point
  vector<lower=-0.5,upper=0.5>[T-1] s_av_rep1; // mean fitness
  
  vector<lower=1,upper=50>[T] c_rep2;               // c at each time point
  vector<lower=2e4,upper=2e9>[T-1] Nb_rep2;       // true bottleneck size at each time point
  vector<lower=-0.5,upper=0.5>[T-1] s_av_rep2; // mean fitness
}
model {
  vector[T] beta_rep1;
  vector[T-1] a_rep1;
  vector[T] beta_rep2;
  vector[T-1] a_rep2;
  real mu;
  real v;
  
  eta_cfus ~ gamma(3.8231428956, 4.1770431245); // from REL606 cfu error parameter fit
  sig ~ pareto(1e-4, 1);

  for (t in 1:T) {
    if (t!=T) {
      s_av_rep1[t] ~ normal(0,0.1);
      a_rep1[t] = exp(-s_av_rep1[t]*6.64/2);
      if (cfus1_rep1[t]!=0){
        cfus1_rep1[t] ~ neg_binomial_2(Nb_rep1[t]*D_rep1[t], Nb_rep1[t]*D_rep1[t]/eta_cfus); // phi = mu/(c-1)
      }
      if (cfus2_rep1[t]!=0){
        cfus2_rep1[t] ~ neg_binomial_2(Nb_rep1[t]*D_rep1[t], Nb_rep1[t]*D_rep1[t]/eta_cfus); // phi = mu/(c-1)
      }
    }
    c_rep1[t] ~ pareto(1, 1);
    beta_rep1[t] = c_rep1[t]/(4*R_rep1[t]);
  }

  for (t in 1:T) {
    if (t!=T) {
      s_av_rep2[t] ~ normal(0,0.1);
      a_rep2[t] = exp(-s_av_rep2[t]*6.64/2);
      if (cfus1_rep2[t]!=0){
        cfus1_rep2[t] ~ neg_binomial_2(Nb_rep2[t]*D_rep2[t], Nb_rep2[t]*D_rep2[t]/eta_cfus); // phi = mu/(c-1)
      }
      if (cfus2_rep2[t]!=0){
        cfus2_rep2[t] ~ neg_binomial_2(Nb_rep2[t]*D_rep2[t], Nb_rep2[t]*D_rep2[t]/eta_cfus); // phi = mu/(c-1)
      }
    }
    c_rep2[t] ~ pareto(1, 1);
    beta_rep2[t] = c_rep2[t]/(4*R_rep2[t]);
  }

  // replicate 1
  for (m in 1:M_rep1) {
    mu = a_rep1[1]*f_rep1[m][1];
    v = (a_rep1[1]*a_rep1[1])*beta_rep1[1] + sig/(4*Nb_rep1[1]*a_rep1[1]*a_rep1[1]);
    for (t in 2:T) {
      f_rep1[m][t] ~ LPTN(mu,sqrt(v + beta_rep1[t]), xi_rep1);
      if (t!=T) {
        mu = a_rep1[t]*(f_rep1[m][t]*v + mu*beta_rep1[t])/(v + beta_rep1[t]);
        v = (a_rep1[t]*a_rep1[t])*(1/((1/v) + (1/beta_rep1[t]))) + sig/(4*Nb_rep1[t]*a_rep1[t]*a_rep1[t]);
      }
    }
  }

  //replicate 2
  for (m in 1:M_rep2) {
    mu = a_rep2[1]*f_rep2[m][1];
    v = (a_rep2[1]*a_rep2[1])*beta_rep2[1] + sig/(4*Nb_rep2[1]*a_rep2[1]*a_rep2[1]);
    for (t in 2:T) {
      f_rep2[m][t] ~ LPTN(mu,sqrt(v + beta_rep2[t]), xi_rep2);
      if (t!=T) {
        mu = a_rep2[t]*(f_rep2[m][t]*v + mu*beta_rep2[t])/(v + beta_rep2[t]);
        v = (a_rep2[t]*a_rep2[t])*(1/((1/v) + (1/beta_rep2[t]))) + sig/(4*Nb_rep2[t]*a_rep2[t]*a_rep2[t]);
      }
    }
  }
}
"""

if len(sys.argv)!=2:
  raise Exception('Specify all files')




expstrain=sys.argv[1].split(',')

strain=expstrain[0]
exp1=expstrain[1]
exp2=expstrain[2][:-1]
days = [0,1,2,3,4]
tv = [str(t) for t in days]

outdir='../finalestimates/'

print('Running HMM in STAN for {}'.format(strain))

#######
freqs_1=pd.read_csv('../{}/data/bc_counts/{}_{}_freqs_filtered_orig.csv'.format(exp1,strain,1))#.set_index('super bc')
Rtot_1=pd.read_csv('../{}/data/bc_counts/{}_{}_Rtot.csv'.format(exp1,strain,1))


cfus_1=pd.read_csv('../{}/{}_cfus.csv'.format(exp1,exp1)).fillna(0)

cond1=cfus_1['Strain']==strain
cond2=cfus_1['Replicate']==1
cfus_1=cfus_1[cond1 & cond2].sort_values(by=['Day'])

fd_1=[]
for i,row in freqs_1.iterrows():
  fd_1.append(list(np.sqrt(row[tv])))
#######
freqs_2=pd.read_csv('../{}/data/bc_counts/{}_{}_freqs_filtered_orig.csv'.format(exp2,strain,2))#.set_index('super bc')
Rtot_2=pd.read_csv('../{}/data/bc_counts/{}_{}_Rtot.csv'.format(exp2,strain,2))


cfus_2=pd.read_csv('../{}/{}_cfus.csv'.format(exp2,exp2)).fillna(0)

cond1=cfus_2['Strain']==strain
cond2=cfus_2['Replicate']==2
cfus_2=cfus_2[cond1 & cond2].sort_values(by=['Day'])

fd_2=[]
for i,row in freqs_2.iterrows():
  fd_2.append(list(np.sqrt(row[tv])))
#######

data_dic ={
  'T': len(days),
  'M_rep1':len(fd_1),
  'xi_rep1':1/(1+np.float(len(fd_1))/500),
  'f_rep1':fd_1,
  'R_rep1':list(np.array(Rtot_1['R'])), #R new
  'D_rep1':list(cfus_1['Dilution']),
	'cfus1_rep1':[np.int(cc) for cc in list(cfus_1['CFU1'])],
  'cfus2_rep1':[np.int(cc) for cc in list(cfus_1['CFU2'])],

  'M_rep2':len(fd_2),
  'xi_rep2':1/(1+np.float(len(fd_2))/500),
  'f_rep2':fd_2,
  'R_rep2':list(np.array(Rtot_2['R'])), #R new
  'D_rep2':list(cfus_2['Dilution']),
  'cfus1_rep2':[np.int(cc) for cc in list(cfus_2['CFU1'])],
  'cfus2_rep2':[np.int(cc) for cc in list(cfus_2['CFU2'])],
}
#print(data_dic)
#print(len(days),len(fd))

sm = pystan.StanModel(model_code=stancode)
fit = sm.sampling(data=data_dic, iter=1000, chains=4)
print(fit)

with open(outdir+'/{}_StanOutput.txt'.format(strain), "w") as text_file:
    print(fit, file=text_file)



sig=pd.DataFrame(fit.extract(['sig'])['sig'],columns=['sig2'])
sig.to_csv(outdir+'/{}_sig.csv'.format(strain),index=False)
eta=pd.DataFrame(fit.extract(['eta_cfus'])['eta_cfus'],columns=['eta'])
eta.to_csv(outdir+'/{}_eta_cfus.csv'.format(strain),index=False)

Nbs=[]
Nes=[]

for rep in [1,2]:
  rept='rep'+str(rep)
  pd.DataFrame(fit.extract(['c_{}'.format(rept)])['c_{}'.format(rept)],columns=tv
      ).to_csv(outdir+'/{}_{}_cs.csv'.format(strain,rep),index=False)


  Nb = pd.DataFrame(fit.extract(['Nb_{}'.format(rept)])['Nb_{}'.format(rept)],columns=tv[:-1])
  har = 1/np.mean(1/np.array(Nb[tv[:-1]]),axis=1)
  arm = np.mean(np.array(Nb[tv[:-1]]),axis=1)
  Nb['harmonic mean'] = har
  Nb['ar mean'] = arm
  Nb.to_csv(outdir+'/{}_{}_Nb.csv'.format(strain,rep),index=False)
  Nbs.append(Nb)

  s = pd.DataFrame(fit.extract(['s_av_{}'.format(rept)])['s_av_{}'.format(rept)],columns=tv[:-1])
  s['av'] = np.mean(np.array(s[tv[:-1]]),axis=1)
  s.to_csv(outdir+'/{}_{}_sav.csv'.format(strain,rep),index=False)

  Ne=Nb.div(np.array(fit.extract(['sig'])['sig']),axis='index')
  Ne.to_csv(outdir+'/{}_{}_Ne.csv'.format(strain,rep),index=False)
  Nes.append(Ne)

pd.DataFrame({
    'harmonic mean':2/(1/np.array(Nbs[0]['harmonic mean']) + 1/np.array(Nbs[1]['harmonic mean'])),
    'ar mean': (np.array(Nbs[0]['ar mean']) + np.array(Nbs[1]['ar mean']))/2,
  }).to_csv(outdir+'/{}_avNb.csv'.format(strain),index=False)


pd.DataFrame({
    'harmonic mean':2/(1/np.array(Nes[0]['harmonic mean']) + 1/np.array(Nes[1]['harmonic mean'])),
    'ar mean': (np.array(Nes[0]['ar mean']) + np.array(Nes[1]['ar mean']))/2,
  }).to_csv(outdir+'/{}_avNe.csv'.format(strain),index=False)
