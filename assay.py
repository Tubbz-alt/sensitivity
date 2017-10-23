import numpy as np

class Assay:
  def __init__(self,params):
    self.params = params

  def throw(self,method):
    #self.parse(method)

    # Gaussian
    if method == 'Gaussian' or 'mu' in self.params.keys():
      if 'limit' in self.params.keys():
        self.params['mu'] = 0
        self.params['sigma'] = self.params['limit']/1.64485
      
      while True:
        value = np.random.normal(self.params['mu'],self.params['sigma'])
        if value >= 0: break

    # Uniform 
    elif method == 'Uniform':
      lowerlimit = 0
      upperlimit = self.params['limit']/0.9
      value = 1.*(upperlimit-lowerlimit)*np.random.rand() + lowerlimit

    # Delta
    elif method == 'Delta':
      value = 1.*self.params['limit']

    # Exponential
    elif method == 'Exponential':
      lam = -math.log(1.-0.9)/self.params['limit']
      value = np.random.exponential(1./lam)

    # Half-cauchy
    elif method == 'Cauchy':
      x0 = 0
      gamma = self.params['limit']/math.tan((0.95-0.5)*math.pi)
      while True:
        value = x0+gamma*np.random.standard_cauchy()
        if value >= 0: break

    # FlatTop
    elif method == 'FlatTop':
      if np.random.rand() > 0.9:
        while(True):
          value = abs(np.random.normal(self.params['limit'], self.params['limit']/1.64))
          #value = abs(np.random.normal(0., self.params['limit']/1.64))
          if value > self.params['limit']: break
      else:
        value = self.params['limit']*np.random.rand()
 
    # Plateau
    elif method == 'Plateau':
      mu = self.params['limit']*np.random.rand()
      sigma = self.params['limit']/1.64
      while True:
        value = np.random.normal(mu,sigma)
        if value >= 0: break

    # Truncated Gaussian
    elif method == 'TruncatedGaussian':
      mu = self.params['original']['mu']
      sigma = self.params['original']['sigma']
      while True:
        value = np.random.normal(mu,sigma)
        if value >= 0: break

    # Truncated Gaussian with Nonzero Mu
    elif method == 'TNZGaussian':
      mu = max(0.,self.params['original']['mu'])
      sigma = self.params['original']['sigma']
      while True:
        value = np.random.normal(mu,sigma)
        if value >= 0: break

    else:
      value = 0

    return value

  def __repr__(self):
    return str(self.params)

from ROOT import * 
import math
gROOT.SetBatch(1)
if __name__ == '__main__':
  assay = Assay({'limit':1.64})
  s = [] 
  histo = TH1D('histo','histo',500,-1,5)
  for i in range(1000000):
    histo.Fill(assay.throw('Plateau'))
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
  


