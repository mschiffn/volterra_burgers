function [pressure_output, N_steps] = burgers_nonlinear_step(pressure_input, f_s, N_samples, delta_z, b, method_interp)
%solve nonlinear part of the Burgers equation
%author: Martin Schiffner
%date: 2009-03-30

tau_z = (1:N_samples) / f_s;        %time axis associated to input signal
N_steps = 0;

while(delta_z > 0)
    %propagation in forward direction
       
    delta_z_act = delta_z;
        
    %check propagation distance delta_z_act to ensure correct ordering
    deriv_est = (pressure_input(2:end) - pressure_input(1:(end - 1))) * f_s;      
    deriv_est_max = max(deriv_est); %maximum of the estimated derivative
        
    step_size_allowed_max = 1 / (b * deriv_est_max);
    
    if delta_z_act >= step_size_allowed_max
            
        %step size is too great, adapt step size
        delta_z_act = step_size_allowed_max / 2;   
    end
    %assertion: delta_z_act < step_size_allowed_max
         
    %commit time-shift for time axis associated to pressure_input
    tau_z_next = tau_z - b * delta_z_act * pressure_input;
        
    %interpolate current waveform linearly to obtain equidistant time axis
    pressure_input = interp1(tau_z_next, pressure_input, tau_z, method_interp, 0);
        
    %reduce propagation distance to go
    delta_z = delta_z - delta_z_act;
    
    %count steps
    N_steps = N_steps + 1;
end
    
while(delta_z < 0)
    %propagation in backward direction
    
    delta_z_act = delta_z;
    
    %check propagation distance delta_z_act to ensure correct ordering
    deriv_est = (pressure_input(2:end) - pressure_input(1:(end - 1))) * f_s;         
    deriv_est_min = min(deriv_est); %minimum of the estimated derivative
    
    step_size_allowed_min = 1 / (b * deriv_est_min);
    
    if delta_z_act <= step_size_allowed_min
            
        %step size is too great in direction of the negative z-axis, adapt step size
        delta_z_act = step_size_allowed_min / 2;   
    end
    %insertion: delta_z_act < step_size_allowed_max
    
    %commit time-shift for time axis associated to pressure_input
    tau_z_next = tau_z - b * delta_z_act * pressure_input;
        
    %interpolate current waveform linearly to obtain equidistant time axis
    pressure_input = interp1(tau_z_next, pressure_input, tau_z, method_interp, 'extrap');
        
    %check for NaN (caused by out of range values) and set them to zero
    pressure_input(isnan(pressure_input)) = 0;
    
    %reduce propagation distance to go
    delta_z = delta_z - delta_z_act;
    
    %count steps
    N_steps = N_steps + 1;
end

%assign changes to output signal
pressure_output = pressure_input;