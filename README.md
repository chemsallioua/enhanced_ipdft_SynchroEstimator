# e-IpDFT SynchroEstimator

The e-IpDFT-SynchroEstimator is a C MEX S-Function based Simulink Block implementation of a Phasor Measurment Unit for Real-Time Simulations that has been developed and tested to comply with industry standards for power measurement. 

## Description
It has been implemented as a low-level C function in Simulink, without reliance on any pre-defined libraries, in order to ensure good performance when running in real-time and compatibility with different Simulink versions. 
The block is customizable and user-friendly, and can be used to measure single-phase or three-phase systems with adjustable parameters such as the acquisition window length, reporting rate, sample rate, and optional timestamping of output data. 

The PMU block uses an iterative e-IpDFT algorithm, which belongs to the DFT (Discrete Fourier Transform) category of algorithms, to estimate synchrophasors (a measure of the phase angle and frequency of an AC voltage or current). This algorithm is known for its balance of accuracy, computational complexity, and response time, and has been extended to compensate for spectral leakage and interharmonic interference.

## Building from source

To compile this S-function, enter
```
mex e_ipDFT_SynchroEstimator.c
```
at the Matlab command line. The mex command compiles and links the e_ipDFT_SynchroEstimator.c file using the default compiler. The mex command creates a dynamically loadable executable for the Simulink software to use.

## Resources
For more info read the following article:
[Open-Source MATLAB-Based PMU Library for HIL Applications Compliant with IEC/IEEE 60255-118-1](https://ieeexplore.ieee.org/document/9978837/)
DOI: [10.1109/AMPS55790.2022.9978837](http://dx.doi.org/10.1109/AMPS55790.2022.9978837)

The Matlab\Simulink Block is available on Mathworks:
[![View Phasor Measurment Unit (Single-Phase M-Class) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/118190-phasor-measurment-unit-single-phase-m-class)
[Single-Phase M-Class PMU (Phasor Measurment Unit)](https://it.mathworks.com/matlabcentral/fileexchange/118190-single-phase-m-class-pmu-phasor-measurment-unit?s_tid=prof_contriblnk)
