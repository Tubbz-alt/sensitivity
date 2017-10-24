# Sensitivity

This is a collection of simple python scripts to calculate the sensitivity for a basic counting experiment under different Bayesian interpretations of radioassay results.

## Purpose

Sensitivity of low background experiments (e.g. double beta decay and dark matter searches) is often calculated by toy Monte Carlo assuming all background components are Poisson random variables. 

The expected background rate (parameter lambda for the Poisson distribution) for each component is essentially the product of the hit efficiency determined by detector simulation and the radioactivity level of the part determined by dedicated radioassay measurements (e.g. by HPGe counting or ICP-MS).

As the radiopurities of the parts improve, the radioactivity levels may not be precisely measured, or worse still, only an upper limit can be set. In such cases, the radioactivity levels are treated as an additional stochastic element in the sensitivity calculation in a Bayesian sense. 

When a radioactivity measurement is made (albeit imprecise), the radioactivity level is typically treated as a Gaussian with mu equals the central value and sigma equals the "error" of the measurement. Whereas only a limit is reported, there isn't an intuitive way to interpret such a result as a Bayesian prior. 

The purpose of this set of scripts is to investigate how the choice of a Bayesian prior for a radioactivity limit affects the calculated sensitivity.

Some ways to interpret a radioactivity limit include:
1. Delta:
2. Gaussian:
  * Allow negative Mu's
  * Mu's are non-negative
3. Uniform:


## Code

### Prerequisites
1. PyROOT (for histogramming and Feldman-Cousins calculator)
2. numpy (for random number generation and miscellaneous numerical utilities)

### Usage
`python toysens.py <ntoys> <true_lambda> <ncomp> <spec_act> <livetime> <method> <nsenstoys>`

1. `ntoys`: Number of toy MCs
2. `true_lambda`: 
3. `ncomp`: 
4. `spec_act`:
5. `livetime`: 
6. `method`: "Truth", "Gaussian",
7. `nsenstoys`:
