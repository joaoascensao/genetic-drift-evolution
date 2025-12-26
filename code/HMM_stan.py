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
  int<lower=0> M;           // # of barcodes
  real<lower=0, upper=1> xi;
  vector[T] f[M]; // sqrt-transformed observed frequencies
  vector[T] R;              // total sequence reads for each time point
  vector[T-1] D;              // cfu dilution rate for each time point, after accounting for the transfer dilution
  int cfus1[T-1];           // CFUs for each time point (not multiplied by dilution rate)
  int cfus2[T-1];           // CFUs for each time point (not multiplied by dilution rate)
}
parameters {
  vector<lower=1,upper=50>[T] c;               // c at each time point
  vector<lower=2e4,upper=2e9>[T-1] Nb;       // true bottleneck size at each time point
  real<lower=1e-4,upper=50> sig;               // sigma^2, descendant number variance
  vector<lower=-0.5,upper=0.5>[T-1] s_av; // mean fitness
  real<lower=1e-4,upper=30> eta_cfus; // error parameter for CFUs - 1
}
model {
  vector[T] beta;
  vector[T-1] a;
  real mu;
  real v;
  
  eta_cfus ~ gamma(3.8231428956, 4.1770431245); // from REL606 cfu error parameter fit

  for (t in 1:T) {
    if (t!=T) {
      s_av[t] ~ normal(0,0.1);
      a[t] = exp(-s_av[t]*6.64/2);
      if (cfus1[t]!=0){
        cfus1[t] ~ neg_binomial_2(Nb[t]*D[t], Nb[t]*D[t]/eta_cfus); // phi = mu/(c-1)
      }
      if (cfus2[t]!=0){
        cfus2[t] ~ neg_binomial_2(Nb[t]*D[t], Nb[t]*D[t]/eta_cfus); // phi = mu/(c-1)
      }
    }
    c[t] ~ pareto(1, 1);
    beta[t] = c[t]/(4*R[t]);
  }

  sig ~ pareto(1e-4, 1);
  for (m in 1:M) {
    mu = a[1]*f[m][1];
    v = (a[1]*a[1])*beta[1] + sig/(4*Nb[1]*a[1]*a[1]);
    for (t in 2:T) {
      f[m][t] ~ LPTN(mu,sqrt(v + beta[t]), xi);
      if (t!=T) {
        mu = a[t]*(f[m][t]*v + mu*beta[t])/(v + beta[t]);
        v = (a[t]*a[t])*(1/((1/v) + (1/beta[t]))) + sig/(4*Nb[t]*a[t]*a[t]);
      }
    }
  }
}
"""

if len(sys.argv)!=6:
  raise Exception('Specify all files')


countdir=sys.argv[1]
meta=pd.read_csv(sys.argv[3])
outdir=sys.argv[2]

exprep=sys.argv[4].split(',')
strain=exprep[0]
replicate=exprep[1][0]
freqs=pd.read_csv(countdir + '/{}_{}_freqs_filtered_orig.csv'.format(strain,replicate))#.set_index('super bc')
Rtot=pd.read_csv(countdir+ '/{}_{}_Rtot.csv'.format(strain,replicate))


cfus=pd.read_csv(sys.argv[5]).fillna(0)

print('Running HMM in STAN for {},{}'.format(strain,replicate))

cond1=meta['Strain']==strain
cond2=meta['Replicate']==np.int(replicate)
mm=meta[cond1 & cond2]


cond1=cfus['Strain']==strain
cond2=cfus['Replicate']==np.int(replicate)
cfus=cfus[cond1 & cond2].sort_values(by=['Day'])

days=np.sort(mm['Day'])
days = days[days>=0]
tv = [str(t) for t in days]

fd=[]
for i,row in freqs.iterrows():
  fd.append(list(np.sqrt(row[tv])))


data_dic ={
  'T': len(days),
  'M':len(fd),
  'xi':1/(1+np.float(len(fd))/500),
  'f':fd,
  'R':list(np.array(Rtot['R'])), #R new
  'D':list(cfus['Dilution']),
	'cfus1':[np.int(cc) for cc in list(cfus['CFU1'])],
  'cfus2':[np.int(cc) for cc in list(cfus['CFU2'])],
}


sm = pystan.StanModel(model_code=stancode)
fit = sm.sampling(data=data_dic, iter=1000, chains=4)
print(fit)

with open(outdir+'/{}_{}_StanOutput.txt'.format(strain,replicate), "w") as text_file:
    print(fit, file=text_file)

pd.DataFrame(fit.extract(['c'])['c'],columns=tv).to_csv(outdir+'/{}_{}_cs.csv'.format(strain,replicate),index=False)

Nb = pd.DataFrame(fit.extract(['Nb'])['Nb'],columns=tv[:-1])
har = 1/np.mean(1/np.array(Nb[tv[:-1]]),axis=1)
arm = np.mean(np.array(Nb[tv[:-1]]),axis=1)
Nb['harmonic mean'] = har
Nb['ar mean'] = arm
Nb.to_csv(outdir+'/{}_{}_Nb.csv'.format(strain,replicate),index=False)

sig=pd.DataFrame(fit.extract(['sig'])['sig'],columns=['sig2'])
sig.to_csv(outdir+'/{}_{}_sig.csv'.format(strain,replicate),index=False)

s = pd.DataFrame(fit.extract(['s_av'])['s_av'],columns=tv[:-1])
s['av'] = np.mean(np.array(s[tv[:-1]]),axis=1)
s.to_csv(outdir+'/{}_{}_sav.csv'.format(strain,replicate),index=False)

Nb.div(np.array(fit.extract(['sig'])['sig']),axis='index').to_csv(outdir+'/{}_{}_Ne.csv'.format(strain,replicate),index=False)


sig=pd.DataFrame(fit.extract(['eta_cfus'])['eta_cfus'],columns=['eta'])
sig.to_csv(outdir+'/{}_{}_eta_cfus.csv'.format(strain,replicate),index=False)