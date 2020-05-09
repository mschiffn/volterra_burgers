function [ pressure_output, N_steps ] = step( pressure_input, f_s, delta_z, a, b, order )
%
% Strang-Marchuk splitting step for the Burgers equation
%
% The Burgers equation is split into a diffusion equation and a nonlinear
% equation.
%
% Split the propagation interval delta_z into two halves (see [1, Sect. 2.2]):
%   1.) apply operator A to the first half
%   2.) apply operator B to thw whole interval
%   3.) apply operator A to the second half
% where the operators A and B equal the diffusion or nonlinear operators.
% 
% INPUT:
%   pressure_input = samples of the acoustic pressure waveform at the boundary
%   f_s = sampling rate (Hz)
%   delta_z = propagation interval
%   a = parameter in Burgers equation
%   b = parameter in Burgers equation
%   order = select the roles of operators A and B
%
% OUTPUT:
%   pressure_output = samples of the acoustic pressure waveform after propagation length delta_z
%   N_steps = 
%
% REFERENCES:
%	[1] I. Faragó, "A modified iterated operator splitting method,"
%       Applied Mathematical Modelling, vol. 32, no. 5B, pp. 1542-1551, Nov. 2008.
%       DOI: 10.1016/j.apm.2007.04.018
%   [2] J. Tavakkoli, D. Cathignol, R. Souchon, O. A. Sapozhnikov, "Modeling of pulsed finite-amplitude focused sound beams in time domain,"
%       J. Acoust. Soc. Am., vol. 104, no. 4, pp. 2061-2072, Oct. 1998
%       DOI: 10.1121/1.423720
%   [3] R. J. Zemp, J. Tavakkoli, R. S. C. Cobbold, "Modeling of nonlinear ultrasound propagation in tissue from array transducers,"
%       J. Acoust. Soc. Am., vol. 113, no. 1, pp. 139-152, Jan. 2003
%       DOI: 10.1121/1.1528926
%
% author: Martin F. Schiffner
% date: 2009-03-30
% modified: 2020-05-08

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 2.) compute splitting step
%--------------------------------------------------------------------------
%
N_samples = numel(pressure_input);
ksi = 0.5;                          %factor used in finite difference scheme (0.5 = Crank-Nicolson method)
method_interp = 'linear';           %method which is used for interpolation in solution to nonlinear problem

%
switch order
    
    case 0
        % solve linear diffusion equation for delta_z / 2 (adjustable finite difference scheme)
        pressure_input = fractional_steps.diffusion( pressure_input, f_s, N_samples, delta_z / 2, a, ksi );
        
        % solve nonlinear equation for delta_z    
        [ pressure_input, N_steps ] = fractional_steps.nonlinear( pressure_input, f_s, N_samples, delta_z, b, method_interp );

        % solve linear diffusion equation for delta_z / 2 (adjustable finite difference scheme)
        pressure_output = fractional_steps.diffusion( pressure_input, f_s, N_samples, delta_z / 2, a, ksi );
     
    case 1
        % solve nonlinear equation for delta_z / 2
        [ pressure_input, N_steps_1 ] = fractional_steps.nonlinear( pressure_input, f_s, N_samples, delta_z / 2, b, method_interp );
        
        % solve linear diffusion equation for delta_z (adjustable finite difference scheme)
        pressure_input = fractional_steps.diffusion( pressure_input, f_s, N_samples, delta_z, a, ksi );
 
        % solve nonlinear equation for delta_z / 2
        [ pressure_output, N_steps_2 ] =  fractional_steps.nonlinear( pressure_input, f_s, N_samples, delta_z / 2, b, method_interp );
        
        N_steps = [N_steps_1, N_steps_2];
        
    otherwise

        % display error message
        error('burgers_frac_steps: invalid argument order in fractional steps scheme');

end

end % function [ pressure_output, N_steps ] = step( pressure_input, f_s, delta_z, a, b, order )
