# Measuring proximity to criticality

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/DLAG/blob/master/LICENSE.md

This project contains a MATLAB implementation of the d<sub>2</sub> measure of proximity to criticality described in Sooter et al., *Defining and measuring proximity to criticality*, bioRxiv (2025). 

<img src="https://github.com/user-attachments/assets/34e8e526-59b4-40e8-8e4a-099933aa4f88" style="width:50%; height:auto;">


Given an observed time series, d<sub>2</sub> quantifies how close the underlying dynamics is to criticality, i.e. how close it is to a scale-invariant boundary between different dynamical regimes. Specifically, d<sub>2</sub> answers the question "How distinguishable (in the information-theoretic sense) is this system from a system at criticality?"

The main code is in the `\src` directory.

The `\demo` directory contains example applications to simulated data and publicly available electrophysiological data from mouse visual cortex [(Senzai et al., Neuron 2019)](https://pubmed.ncbi.nlm.nih.gov/30635232/).

The software has been tested with MATLAB version R2024b.

## Citing this work

Please cite [this reference] and [this codepack] if you use any portion of this code in your own work.

## Installation guide

Simply download and extract the latest release of this project to your local working directory. You will also need to install the 
[Optimization Toolbox](https://www.mathworks.com/help/optim/index.html) and [Parallel Computing Toolbox](https://www.mathworks.com/help/parallel-computing/index.html) if these are not already included in your Matlab build.

## Contact
For questions, please contact Sam Sooter at jssooter@gmail.com or Woodrow Shew at shew@uark.edu

## License
[MIT](LICENSE.md)
