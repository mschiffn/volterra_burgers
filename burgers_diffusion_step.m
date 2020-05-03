function pressure_output = burgers_diffusion_step(pressure_input, f_s, N_samples, delta_z, a, ksi)
%finite difference scheme for Burgers equation
%author: Martin Schiffner
%date: 2009-03-30

%properties for delta_z > 0:
%ksi = 0: explicit method           ,global error = O(delta_z) + O(1/f_s^2)     ,stability must be verified (r <= 0.5)
%ksi = 0.5: Crank-Nicolson method   ,global error = O(delta_z^2) + O(1/f_s^2)   ,always stable (ksi >= 0.5)
%ksi = 1: implicit method           ,global error =   ?                         ,always stable (ksi >= 0.5)

%properties for delta_z < 0:    ksi_prime = 1 - ksi
%ksi = 1: explicit method           ,never stable for r > 0
%ksi = 0.5: Crank-Nicolson method   ,never stable for r > 0
%ksi = 0: implicit method           ,never stable for r > 0

r = a * abs(delta_z) * f_s^2;       %compute modulus (if ksi < 0.5: stability check is necessary)

J = diag(2 * ones(1, N_samples)) - diag(ones(1, N_samples - 1), 1) - diag(ones(1, N_samples - 1), -1);
B = eye(N_samples) + r * ksi * J;
C = eye(N_samples) - r * (1 - ksi) * J;

%check direction of propagation
if delta_z > 0
    
    %propagation in forward direction   
    
    %solve linear equations (assumption: pressure vanishes on the borders)
    pressure_output = (B \ (C * pressure_input'))';

elseif delta_z < 0
    
    %propagation in backward direction
    
    %solve linear equations (assumption: pressure vanishes on the borders)
    pressure_output = (C \ (B * pressure_input'))';
end
