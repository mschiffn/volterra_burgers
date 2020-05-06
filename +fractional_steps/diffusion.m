function pressure_output = diffusion( pressure_input, f_s, N_samples, delta_z, a, xi )
%
% finite difference time domain (FDTD) method for
% the diffusion equation
%
%--------------------------------------------------------------------------
% properties for delta_z > 0
%--------------------------------------------------------------------------
% xi        method          global error                    stability
% 0         explicit        O(delta_z) + O(1/f_s^2)         must be verified (r <= 0.5)
% 0.5       Crank-Nicolson	O(delta_z^2) + O(1/f_s^2)       always stable (xi >= 0.5)
% 1         implicit        ?                               always stable (xi >= 0.5)
%
%--------------------------------------------------------------------------
% properties for delta_z < 0:    ksi_prime = 1 - xi
%--------------------------------------------------------------------------
% xi        method          global error                    stability
% 1         explicit                                        never stable for r > 0
% 0.5       Crank-Nicolson                                  never stable for r > 0
% 0         implicit                                        never stable for r > 0
%
% author: Martin Schiffner
% date: 2009-03-30
% modified: 2020-05-06

	%----------------------------------------------------------------------
	% 1.) check arguments
	%----------------------------------------------------------------------
	% ensure valid number of input arguments
	narginchk( 6, 6 );

    %----------------------------------------------------------------------
    % 2.) FDTD method
    %----------------------------------------------------------------------
    % specify matrices
    r = a * abs( delta_z ) * f_s^2;	% compute modulus (if xi < 0.5: stability check is necessary)

    J = 2 * eye( N_samples ) - diag( ones( 1, N_samples - 1 ), 1 ) - diag( ones( 1, N_samples - 1 ), -1 );
    B = eye( N_samples ) + r * xi * J;
    C = eye( N_samples ) - r * ( 1 - xi ) * J;

	% check direction of propagation
	if delta_z > 0

        %------------------------------------------------------------------
        % a) forward propagation
        %------------------------------------------------------------------
        % solve linear system (assumption: pressure vanishes on the borders)
        pressure_output = ( B \ ( C * pressure_input' ) )';

	elseif delta_z < 0

        %------------------------------------------------------------------
        % b) backward propagation
        %------------------------------------------------------------------
        % solve linear system (assumption: pressure vanishes on the borders)
        pressure_output = ( C \ ( B * pressure_input' ) )';

    end % if delta_z > 0

end % function pressure_output = diffusion( pressure_input, f_s, N_samples, delta_z, a, xi )
