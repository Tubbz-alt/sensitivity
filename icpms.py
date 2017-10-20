# ICPMS simulator
# 
# heiman.tsang@pnnl.gov

import os
import cPickle as pickle
import numpy as np
from ROOT import *
from TimeProfiler import *

gROOT.SetBatch(1)
#gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)

# Abstraction of an ICPMS machine: Converts concentration to number of counts following Poisson statistics.
class ICPMSMachine:
  def __init__(self,prob):
    self.prob = float(prob)   # counts/ppt  (200-300 counts/(pg/V))
  def __repr__(self):
    return str({'param':self.prob})
  def measure(self,sample,num):
    d = {}
    for name in sample.impurities.keys():
      d[name] = [np.random.poisson(self.prob*sample.impurities[name]) for i in range(num)]
    return d

# Bulk material from a vendor with normal-distributed impurity levels.
class Material:
  def __init__(self,name,impurities):
    self.name = name
    self.impurities = impurities  # {"name": [mean, sigma]}  in pg/g for samples, pg for blanks and tracers
  def __repr__(self):
    return str([self.name, self.impurities])
  def create_sample(self,mass):
    d = {}
    for name in self.impurities.keys():
      while 1:
        d[name] = (mass if mass > 0 else 1.)*np.random.normal(self.impurities[name][0],self.impurities[name][1])
        if d[name] >= 0: break
    return Sample(mass,d)

# Sample taken from a "Material"
class Sample:
  def __init__(self,mass,impurities):
    self.mass = mass # dissolved mass in g. Set to 0 for blanks and tracers.
    self.impurities = impurities  # {"name": value in pg}

  def __repr__(self):
    return str({'mass': self.mass, 'impurities':self.impurities})

  def __add__(self,that):
    x = self.impurities
    y = that.impurities
    return Sample(self.mass+that.mass, { k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y) })

# Simulates the assay procedure from sample taking to analysis result.
#@profile
def assay(material, blank_mat, tracer_mat, **options):

  n_samples = options['n_samples'] if 'n_samples' in options else 3
  n_blanks = options['n_blanks'] if 'n_blanks' in options else 3
  icpms_prob = options['icpms_prob'] if 'icpms_prob' in options else 200.
  n_repetitions = options['n_repetitions'] if 'n_repetitions' in options else 3
  mass = options['mass'] if 'mass' in options else 1.

  # Create samples
  samples = []
  blanks = []

  for i in range(n_samples):
    a = material.create_sample(mass) #0.05)
    b = blank_mat.create_sample(0)
    c = tracer_mat.create_sample(0)
    samples.append(a+b+c)

  for i in range(n_blanks):
    b = blank_mat.create_sample(0)
    c = tracer_mat.create_sample(0)
    blanks.append(b+c)

  # Measure
  machine = ICPMSMachine(icpms_prob)

  sample_meas = [machine.measure(s,n_repetitions) for s in samples]
  blank_meas  = [machine.measure(b,n_repetitions) for b in blanks]

  sample_masses = [s.mass for s in samples]

  if 'verbose' in options and options['verbose'] > 0:
    print('ICPMS:', machine)
    print('Samples:')
    for s in samples: print(s)
    print('Blanks:')
    for s in blanks: print(s)
    print('Sample measurement:')
    for s in sample_meas: print(s)
    print('Blank measurement:')
    for s in blank_meas: print(s)

  return analyze(sample_meas, blank_meas, tracer_mat, sample_masses)

# ICPMS data analysis. 
# Input: count rates of samples and blanks from ICPMS, tracer amount.
# Output: impurity concentrations.
#@profile
def analyze(samples, blanks, tracer_mat, masses):
  #print(samples)
  #print(blanks)
  #print(tracer_mat)
  
  tracer_for = {'U238':'U233'} #,'Th232':'Th229'}

  ans = {}
  fullresult = {}
  lod = {}
  for species in tracer_for.keys():
    tracer = tracer_for[species]
    # y-hat
    species_est = [1.* sum(s[species]) / sum(s[tracer]) * tracer_mat.impurities[tracer][0] for s in samples]
    # b-hat 
    blank_est   = [1.* sum(b[species]) / sum(b[tracer]) * tracer_mat.impurities[tracer][0] for b in blanks] 
    mean_blank_est = np.mean(blank_est)
    #print(species, blank_est)
    #print(species, species_est)
    # x-hat 
    net_species_est = [y - mean_blank_est for y in species_est]
    # normalize for dissolved mass
    net_species_est = [net_species_est[i]/masses[i] for i in range(len(net_species_est))]
    # mu_s-hat
    mean_net_species_est = np.mean(net_species_est)
    
    ans[species] = mean_net_species_est
    fullresult[species] = list(net_species_est)
    # assert same dissolved masses
    lod[species] = 3.*np.std(blank_est)/np.mean(masses)

  return ans, fullresult, lod

@profile
def toy(mu_s, sigma_s, mu_b, sigma_b, n_iter):

  copper_mat = Material('copper',{'U238':[mu_s, mu_s*sigma_s]}) #,'Th232':[5,1]})
  blank_mat  = Material('blank', {'U238':[mu_b, mu_b*sigma_b]}) #,'Th232':[5,1]})
  tracer_mat = Material('tracer',{'U233':[0.5,0.]}) #, 'Th229':[1,0]})

  ans = []
  for i in range(n_iter):
    #if i % (n_iter // 10) == 0: print(i)
    result, fullresult, lod = assay(copper_mat, blank_mat, tracer_mat)#, verbose=1 if i == 0 else 0)
    ans.append({'result': result['U238'], 'fullresult': fullresult['U238'], 'lod': lod['U238']})

  metadata = {
    'mu_s': mu_s,
    'sigma_s': sigma_s,
    'mu_b': mu_b,
    'sigma_b': sigma_b,
    'niter': n_iter,
  }

  data = {
    'metadata': metadata,
    'data': ans,
  }

  return data

# Main program
@profile
def main():
  #mus_range =    [1.] 
  #mub_range =    [10.e-2]
  #sigmas_range = [10.e-2]
  #sigmab_range = [10.e-2]
  mus_range =    [0., 1.e-2, 2.e-2, 5.e-2, 10.e-2, 20.e-2, 50.e-2, 1.]
  mub_range =    [0., 1.e-2, 2.e-2, 5.e-2, 10.e-2, 20.e-2, 50.e-2, 1.]
  sigmas_range = [0., 1.e-2, 2.e-2, 5.e-2, 10.e-2, 20.e-2]
  sigmab_range = [0., 1.e-2, 2.e-2, 5.e-2, 10.e-2, 20.e-2]
  niter = 100000

  for isigmas,sigmas in enumerate(sigmas_range):
    for isigmab,sigmab in enumerate(sigmab_range):

      data_folder = 'data2/sigma_%.3f_%.3f' % (sigmas, sigmab)
      os.system('mkdir -p %s' % data_folder)
  
      for imus,mus in enumerate(mus_range):
        for imub,mub in enumerate(mub_range):
    
              print(mus,sigmas,mub,sigmab,niter)
              fn = '%s/%.3f_%.3f_%.3f_%.3f.pickle' % (data_folder,mus,sigmas,mub,sigmab)
    
              if os.path.isfile(fn):
                print('File exists. Skip.')
                continue
    
              data = toy(mus,sigmas,mub,sigmab,niter)
              pickle.dump(data,open(fn,'wb'))
              print_prof_data()
              clear_prof_data()
    
if __name__ == '__main__':
  main()


