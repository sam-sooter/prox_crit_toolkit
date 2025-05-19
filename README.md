# Measuring proximity to criticality

MATLAB and Python implementation of the d<sub>2</sub> measure of proximity to criticality described in Sooter et al. [link to paper]

<img src="https://github.com/user-attachments/assets/d1ef48e4-b77b-49a2-9905-53daa52da9cb" style="width:50%; height:auto;">

Given an observed time series, d<sub>2</sub> quantifies how close the underlying dynamics is to criticality, i.e. how close it is to a scale-invariant boundary between different dynamical regimes. Specifically, d<sub>2</sub> answers the question "How distinguishable (in the information-theoretic sense) is this system from a system at criticality?" This approach goes beyond previous methods for assessing proximity to criticality in 

The `\src\matlab` and `\src\python` directories contain the MATLAB and python implementations of d<sub>2</sub>, respectively.

The `\test` directory contains example applications of our approach on simulated data and electrophysiological data from mouse visual cortex [Buzsaki].

The software has been tested with MATLAB version xx-yy and Python version aa-bb. 
## Citing this work

For proper attribution, please cite this reference and [this codepack] if
you use any portion of this code in your own work.

## System requirements


## Installation guide

You may need to install the 
[Matlab Optimization Toolbox]
before getting started, if it's not already installed in your Matlab build.

[add...]

## Getting started

Execute `startup.m` to add all necessary dependencies to the Matlab path.
Then to get familiar with the methods in this codepack and their usage, check out the
[`demo`](demo) directory. There, ??

## Project wiki

For additional information on getting started, as well as subtler usage details, see
this project's wiki [link]

## Contact
For questions, please contact Sam Sooter at sooter@uark.edu or Woodrow Shew at shew@uark.edu

## License
[MIT](LICENSE.md)
