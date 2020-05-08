# Evaluation of an Analytical Volterra Series Solution to the Burgers Equation

![GitHub repo size][size-image]
![GitHub All Releases][downloads-image]

[size-image]: https://img.shields.io/github/repo-size/mschiffn/volterra_burgers
[downloads-image]: https://img.shields.io/github/downloads/mschiffn/volterra_burgers/total

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
both methods by evaluating
various error metrics.
It additionally creates
a movie illustrating
the deformation of
the wave and
the accompanying generation of
harmonics.

The package +volterra contains
the functions for
the proposed Volterra polynomial, whereas
the package +fractional_steps contains
the functions for
the fractional steps reference method.

## Results

The following animation depicts
the result of
the 10th-degree Volterra polynomial.

![Animation](./burgers_propagation_hp.gif)

## References :notebook:

The solution and exemplary images were published in:

1. M. Schiffner, M. Mleczko, and G. Schmitz, "Evaluation of an Analytical Solution to the Burgers Equation Based on Volterra Series", 2009 IEEE Int. Ultrasonics Symp. (IUS), pp. 573--576, [![DOI:10.1109/ULTSYM.2009.0096](https://img.shields.io/badge/DOI-10.1109%2FULTSYM.2009.0096-blue)](http://dx.doi.org/10.1109/ULTSYM.2009.0096)
2. M. Schiffner, M. Mleczko, and G. Schmitz, "Application of Volterra Series to Ultrasound Imaging", NAG/DAGA 2009 Int. Conf. Acoustics, Rotterdam, Mar. 2009, pp. 301--304
3. M. Schiffner, M. Mleczko, and G. Schmitz, "Application of Volterra Series to the Detection of Ultrasound Contrast Agents", World Congr. Medical Physics and Biomedical Engineering, Sep. 7--12, 2009, Munich, Germany. IFMBE Proceedings, vol. 25/2, pp. 478--481 [![DOI:10.1007/978-3-642-03879-2_134](https://img.shields.io/badge/DOI-10.1007%2F978--3--642--03879--2__134-blue)](http://dx.doi.org/10.1007/978-3-642-03879-2_134)
