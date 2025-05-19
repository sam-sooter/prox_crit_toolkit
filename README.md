# Measuring proximity to criticality

MATLAB and Python implementation of the d<sub>2</sub> measure of proximity to criticality described in Sooter et al. [link to paper]

<img src="https://github.com/user-attachments/assets/d1ef48e4-b77b-49a2-9905-53daa52da9cb" style="width:50%; height:auto;">

Given an observed time series, d<sub>2</sub> quantifies how close the underlying dynamics is to criticality, i.e. how close it is to a scale-invariant boundary between different dynamical regimes. Specifically, d<sub>2</sub> answers the question "How distinguishable (in the information-theoretic sense) is this system from a system at criticality?" This approach goes beyond previous methods for assessing proximity to criticality in 

The \src\matlab directory contains the  to compute d<sub>2</sub>, with MATLAB code in the \src\matlab sudirector

The \test directory contains example applications of our approach on simulated data and electrophysiological data from mouse visual cortex [Buzsaki].

The software has been tested with MATLAB version xx-yy and Python version aa-bb. The MATLAB implementation requires the Optimization Toolbox and [...].
