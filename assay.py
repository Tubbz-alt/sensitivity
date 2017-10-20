import numpy as np

class Assay:
  def __init__(self,params):
    self.params = params

  def parse(self,method):
    # Gaussian
    if method == 'Gaussian' or 'mu' in self.params.keys():
      if 'limit' in self.params.keys():
        self.params['mu'] = 0
        self.params['sigma'] = self.params['limit']/1.64485

    # Uniform 
    elif method == 'Uniform':
      self.params['lowerlimit'] = 0
      self.params['upperlimit'] = self.params['limit']/0.9

    # Delta
    #elif method == 'Delta':

    # Exponential
    elif method == 'Exponential':
      self.params['lambda'] = -math.log(1.-0.9)/self.params['limit']

    # Half-cauchy
    elif method == 'Cauchy':
      self.params['x0'] = 0
      self.params['gamma'] = self.params['limit']/math.tan((0.95-0.5)*math.pi)

    # FlatTop
    #elif method == 'FlatTop':
      

  def throw(self,method):
    self.parse(method)

    # Gaussian or Half-gaussian
    if method == 'Gaussian' or 'mu' in self.params.keys():
      
      while True:
        value = np.random.normal(self.params['mu'],self.params['sigma'])
        if value >= 0: break

    # Uniform 
    elif method == 'Uniform':
      value = 1.*(self.params['upperlimit']-self.params['lowerlimit'])*np.random.rand() + self.params['lowerlimit']

    # Delta
    elif method == 'Delta':
      value = self.params['limit']

    # Exponential
    elif method == 'Exponential':
      value = np.random.exponential(1./self.params['lambda'])

    # Half-cauchy
    elif method == 'Cauchy':
      while True:
        value = self.params['x0']+self.params['gamma']*np.random.standard_cauchy()
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
    histo.Fill(assay.throw('FlatTop'))
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
  


