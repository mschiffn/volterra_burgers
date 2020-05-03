function [pressure_output, N_steps] = burgers_frac_steps_m(pressure_input, f_s, delta_z, a, b, order)
%author: Martin Schiffner
%date: 2009-03-30

N_samples = numel(pressure_input);
ksi = 0.5;                          %factor used in finite difference scheme (0.5 = Crank-Nicolson method)
method_interp = 'linear';           %method which is used for interpolation in solution to nonlinear problem

switch order
    
    case 0
        % solve linear diffusion equation for delta_z / 2 (adjustable finite difference scheme)
        pressure_input = burgers_diffusion_step(pressure_input, f_s, N_samples, delta_z / 2, a, ksi);
        
        % solve nonlinear equation for delta_z    
        [pressure_input, N_steps] = burgers_nonlinear_step(pressure_input, f_s, N_samples, delta_z, b, method_interp);

        % solve linear diffusion equation for delta_z / 2 (adjustable finite difference scheme)
        pressure_output = burgers_diffusion_step(pressure_input, f_s, N_samples, delta_z / 2, a, ksi);
     
    case 1
        % solve nonlinear equation for delta_z / 2
        [pressure_input, N_steps_1] = burgers_nonlinear_step(pressure_input, f_s, N_samples, delta_z / 2, b, method_interp);
        
        % solve linear diffusion equation for delta_z (adjustable finite difference scheme)
        pressure_input = burgers_diffusion_step(pressure_input, f_s, N_samples, delta_z, a, ksi);
 
        % solve nonlinear equation for delta_z / 2
        [pressure_output, N_steps_2] =  burgers_nonlinear_step(pressure_input, f_s, N_samples, delta_z / 2, b, method_interp);
        
        N_steps = [N_steps_1, N_steps_2];
        
    otherwise
        % display error message
        error('burgers_frac_steps: invalid argument order in fractional steps scheme');     
end
