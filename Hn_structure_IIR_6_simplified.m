function y_out = Hn_structure_IIR_6_simplified(input, filter_fft, filter_fft_2, f_s, a, b, y_out_prec)
%Burgers equation:
%Use simplified structure of n-th-order Volterra kernel with arbitrary
%linear FIR filters to compute the output signal in time-domain
%input: arbitrary discrete-time input signal
%f: impulse response of second FIR filter to use in the structure
%   (exp(a*s^2*x))
%f_s: sampling rate
%a,b: the fluid's parameters
%note: the same sampling rate f_s is assumed for input and f
%
%author: Martin Schiffner
%date: 2008-09-28

%compute output of n-th-order homogeneous Volterra-System recursively
y_out = Hn_structure_IIR_recursion_6_simplified({input, input}, filter_fft, filter_fft_2, f_s, a, b, y_out_prec);
