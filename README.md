# Evaluation of an Analytical Volterra Series Solution to the Burgers Equation

![GitHub repo size][size-image]
![GitHub All Releases][downloads-image]
[![View on File Exchange][matlab-fex-image]][matlab-fex-url]

[size-image]: https://img.shields.io/github/repo-size/mschiffn/volterra_burgers
[downloads-image]: https://img.shields.io/github/downloads/mschiffn/volterra_burgers/total
[matlab-fex-image]: https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg
[matlab-fex-url]: https://de.mathworks.com/matlabcentral/fileexchange/75413-volterra-series-for-the-burgers-equation

Comparison of
a Volterra series solution to
the Burgers equation to
a fractional steps method.

## Background

The Burgers equation models
the nonlinear propagation of
plane ultrasonic waves in
homogeneous viscous fluids.

Its incorporation into
fast tissue harmonic imaging or
the detection of
ultrasound contrast agents potentially improves
these imaging modes.
Moreover,
the decompositions of
arbitrary types of
waves into
steered plane waves permit
the application of
this model
the Burgers equation to
other types of waves.

## Content

The script "eval_methods.m" compares
both methods for
a simple example by evaluating
various error metrics.
It additionally creates
a short movie illustrating
the wave propagation.

The package +volterra contains
the functions for
the proposed Volterra polynomial, whereas
the package +fractional_steps contains
the functions for
the fractional steps reference method.

## Results

The following animation depicts
the waveforms predicted by
the 10th-degree Volterra polynomial for
distilled water
(propagation distance: 25 cm, step length: 0.5 cm).
The initial waveform was
a modulated Gaussian pulse
(center frequency: 3.5 MHz, amplitude: 750 kPa).
The deformation of
the shape (see left plot) involved
the generation of
harmonics (see right plot), which were subject to
significantly increased absorption.

![Propagation of modulated Gaussian pulse](./burgers_propagation_hp.gif)

## References :notebook:

The solution and exemplary applications were published in:

1. M. Schiffner, M. Mleczko, and G. Schmitz, "Evaluation of an analytical solution to the Burgers equation based on Volterra series," in 2009 IEEE Int. Ultrasonics Symp. (IUS), Rome, Sep. 2009, pp. 573-576, [![DOI:10.1109/ULTSYM.2009.5442057](https://img.shields.io/badge/DOI-10.1109%2FULTSYM.2009.5442057-blue)](http://dx.doi.org/10.1109/ULTSYM.2009.5442057)
2. M. Schiffner, M. Mleczko, and G. Schmitz, "Application of Volterra Series to Ultrasound Imaging," in NAG/DAGA 2009 Int. Conf. Acoustics, Rotterdam, Mar. 2009, pp. 301--304
3. M. Schiffner, M. Mleczko, and G. Schmitz, "Application of Volterra Series to the Detection of Ultrasound Contrast Agents," in World Congr. Medical Physics and Biomedical Engineering, Sep. 2009, Munich. IFMBE Proceedings, vol. 25/2, pp. 478-481, [![DOI:10.1007/978-3-642-03879-2_134](https://img.shields.io/badge/DOI-10.1007%2F978--3--642--03879--2__134-blue)](http://dx.doi.org/10.1007/978-3-642-03879-2_134)
