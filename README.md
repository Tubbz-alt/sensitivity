# Sensitivity

This is a collection of simple python scripts to calculate the sensitivity for a basic counting experiment under different Bayesian interpretations of radioassay results.

## Purpose

Sensitivity of low background experiments (e.g. double beta decay and dark matter searches) is often calculated by toy Monte Carlo assuming all background components are Poisson random variables. 

The expected background rate (parameter lambda for the Poisson distribution) for each component is essentially the product of the hit efficiency determined by detector simulation and the radioactivity level of the part determined by dedicated radioassay measurements (e.g. by HPGe counting or ICP-MS).

As the radiopurities of the parts improve, the radioactivity levels may not be precisely measured, or worse still, only an upper limit can be set. In such cases, the radioactivity levels are treated as an additional stochastic element in the sensitivity calculation in a Bayesian sense. 

When a radioactivity measurement is made (albeit imprecise), the radioactivity level is typically treated as a Gaussian with mu equals the central value and sigma equals the "error" of the measurement. Whereas only a limit is reported, there isn't an intuitive way to interpret such a result as a Bayesian prior. 

The purpose of this set of scripts is to investigate how the choice of a Bayesian prior for a radioactivity limit affects the calculated sensitivity, and hence inform us the proper choice.

## Details

The user code is in `toysens.main()`.

First, the user needs to specify, for each part of the Detector in the actual experiment, 
the true specific activity, the hit efficiency and the mass. For example, for `ncomp` identical parts:
```
from detector import Detector, Component
comps = []
for ic in range(ncomp):
  comps.append(Component(spec_act,eff,mass))
det = Detector(comps)
```
Then define the detector for radioassay, e.g. an HPGe detector.:
```
from facilities import HPGe
ge = HPGe(10./86400.)
```
Let the detector perform an assay measurement and associate it with the assayed component:
```
from assay import Assay
result = ge.count(comp.trueimp,1,livetime=14*86400)
comp.assay = Assay(result)
```
Finally perform a sensitivity calculation:
```
livetime = 10*365*86400 # Experiment livetime in seconds
method = 'Uniform'      # Way to interpret limits
nsenstoys = 10000       # Number of toys used in sensitivity calculation
usetruth = False        # Use true radioassay values. If False, use assayed values and interpret according to "method".
sens, uls = calc_sens(det, method, livetime, nsenstoys, usetruth)
```
It will return the 90% C.L. upper limits (`uls`) and their median (`sens`)

Some ways to interpret a radioactivity limit include:
1. Delta:
2. Gaussian:
   1. Allow negative Mu's
   2. Mu's are non-negative
3. Uniform:

## Code

### Prerequisites
1. Python 2.7
2. PyROOT (for histogramming and Feldman-Cousins calculator)
3. numpy (for random number generation and miscellaneous numerical utilities)

### Usage
`python toysens.py <ntoys> <true_lambda> <ncomp> <spec_act> <livetime> <method> <nsenstoys>`

1. `ntoys`: Number of toy MCs
2. `true_lambda`: 
3. `ncomp`: 
4. `spec_act`:
5. `livetime`: 
6. `method`: "Truth", "Gaussian",
7. `nsenstoys`:
