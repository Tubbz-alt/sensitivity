import numpy as np

class Assay:
  def __init__(self,params):
    self.params = params

  def throw(self,method):
    #self.parse(method)
    #print('assay.throw method',method)

    if 'mu' in self.params.keys():  # Discovery
      if method in ['Central', 'DeltaUniform', 'DeltaGauss0', 'CTG']:  # Delta/*
        value = 1.*self.params['mu']
      elif method in ['Delta', 'Uniform', 'Gaussian', 'TruncatedGaussian']:  # Gauss/*
        while True:
          value = np.random.normal(self.params['mu'],self.params['sigma'])
          if value >= 0: break
      elif method[:2] == 'V:': 
        nsig = float(method.split(':')[1])
        nparts = float(method.split(':')[2])
        value = self.params['original']['mu'] + nsig*self.params['original']['sigma']/math.sqrt(nparts)
        value = max(value,0)
      elif method == 'FullGauss': 
        value = np.random.normal(self.params['mu'],self.params['sigma'])
      else:
        value = 0

    else:  # Limit

      if method in ['Central', 'Delta']:  #*/Delta
        value = 1.*self.params['limit']

      elif method in ['DeltaUniform', 'Uniform']: #*/Uniform
        lowerlimit = 0
        upperlimit = self.params['limit']/0.9
        value = 1.*(upperlimit-lowerlimit)*np.random.rand() + lowerlimit

      elif method in ['DeltaGauss0', 'Gaussian']: #*/Gauss0
        mu = 0
        sigma = self.params['limit']/1.64485
        while True:
          value = np.random.normal(mu,sigma)
          if value >= 0: break

      elif method in ['CTG','TruncatedGaussian']:  #*/Gauss
        mu = self.params['original']['mu']
        sigma = self.params['original']['sigma']
        while True:
          value = np.random.normal(mu,sigma)
          if value >= 0: break

      elif method[:2] == 'V:':
        nsig = float(method.split(':')[1])
        nparts = float(method.split(':')[2])
        value = self.params['original']['mu'] + nsig*self.params['original']['sigma']/math.sqrt(nparts)
        value = max(value,0)

      elif method == 'FullGauss': 
        mu = self.params['original']['mu']
        sigma = self.params['original']['sigma']
        value = np.random.normal(mu,sigma)

      else:
        value = 0

    #print('assay',value)
    return value

  def __repr__(self):
    return str(self.params)

from ROOT import * 
import math
from facilities import HPGe
gROOT.SetBatch(1)
if __name__ == '__main__':
  ge = HPGe(10./86400.)
  
  assay = Assay(ge.count(100e-6,1,livetime=14*86400))
 
  s = [] 
  histo = TH1D('histo','histo',500,-10e-6,200e-6)
  for i in range(100000):
    histo.Fill(assay.throw('Gaussian'))
    #histo.Fill(assay.throw('Gaussian'))
    #histo.Fill(assay.throw('Delta'))
    #histo.Fill(assay.throw('Uniform'))

  #histo.Scale(1./0.01)

  canvas = TCanvas('canvas','canvas',1024,768)
  canvas.SetGridx()
  canvas.SetGridy()
  '''
  fGaus = TF1('fGaus','[0]/sqrt(2*pi)/[2] * exp(-((x-[1])**2)/[2]/[2]/2)',0,2) 
  #fGaus = TF1('fGaus','gaus',0,2) 
  fGaus.SetParameter(0,1)
  fGaus.SetParameter(1,0)
  fGaus.SetParameter(2,1)
  histo.Fit(fGaus,'r')
  '''
  histo.Draw()
  canvas.SaveAs('test.png')
  


