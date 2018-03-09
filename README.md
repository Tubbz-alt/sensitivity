# Sensitivity

This is a collection of simple python scripts to calculate the sensitivity for a basic counting experiment under different interpretations of radioassay results.

## Purpose

Sensitivity of low background experiments (e.g. double beta decay and dark matter searches) is often calculated by toy Monte Carlo assuming all background components are Poisson random variables. 

The expected background rate (parameter lambda for the Poisson distribution) for each component is essentially the product of the hit efficiency determined by detector simulation and the radioactivity level of the part determined by dedicated radioassay measurements (e.g. by HPGe counting or ICP-MS).

As the radiopurities of the parts improve, the radioactivity levels may not be precisely measured, sometime only an upper limit can be set. In such cases, the radioactivity levels are treated as a Bayesian prior. 

The purpose of this set of scripts is to investigate how the choice of a Bayesian prior for a radioactivity limit affects the calculated sensitivity, and hence inform us the proper choice.

## Code

### Prerequisites
1. Python 2.7
2. PyROOT (for histogramming and Feldman-Cousins calculator)
3. numpy (for random number generation and miscellaneous numerical utilities)

### Usage

The user code is in `toysens.main()`.

First, the user needs to specify, for each part of the Detector in the actual experiment, 
the true specific activity, the hit efficiency and the mass. For example, for `ncomp` identical parts:
```
from detector import Detector, Component
comps = []
spec_act = 1e-6 # Specific activity of the part in Bq/kg
eff = 0.1       # Hit efficiency of the part in the actual detector
mass = 1.       # Mass of the part in kg
for ic in range(ncomp):
  comps.append(Component(spec_act,eff,mass))
det = Detector(comps)
```
Then define the detector for radioassay, e.g. an HPGe detector.:
```
from facilities import HPGe
ge = HPGe(10./86400.)  # Background rate of HPGe in counts per second
```
Let the detector perform an assay measurement and associate it with the assayed component:
```
from assay import Assay
assayed_mass = 1.   # Mass of part being assayed in kg
result = ge.count(comp.trueimp,assayed_mass,livetime=14*86400)
comp.assay = Assay(result)
```
Here `Assay` performs random draws from the selected Bayesian prior.

Finally perform a sensitivity calculation:
```
livetime = 10*365*86400 # Experiment livetime in seconds
method = 'Uniform'      # Way to interpret limits
nsenstoys = 10000       # Number of toys used in sensitivity calculation
usetruth = False        # If True, use true radioassay values. If False, use assayed values and interpret according to "method".
sens, uls = calc_sens(det, method, livetime, nsenstoys, usetruth)
```
It will return the 90% C.L. upper limits (`uls`) and their mean (`sens`)

Available priors (`method`): 

| Option | Central value | Limits |
| ------ | ------------- | ------ |
| `Central` | Delta | Delta |
| `DeltaUniform` | Delta | Uniform |
| `DeltaGauss0` | Delta | Gaussian at 0 |
| `CTG` | Delta | Truncated Gaussian |
| `Delta` | Gaussian | Delta |
| `Uniform` | Gaussian | Uniform |
| `Gaussian` | Gaussian | Gaussian at 0 |
| `TruncatedGaussian` | Gaussian | Truncated Gaussian |

(See assay.py)

### Quick start

`python toysens.py <ntoys> <true_lambda> <ncomp> <spec_act> <livetime> <method> <nsenstoys>`

1. `ntoys`: Number of toy MCs for the assay procedure
2. `true_lambda`: Expected true total background rate
3. `ncomp`: Number of parts
4. `spec_act`: Specific activity of every part
5. `livetime`: Livetime of the experiment
6. `method`: Method to interpret assay results
7. `nsenstoys`: Number of toy MCs used in each sensitivity calculation

