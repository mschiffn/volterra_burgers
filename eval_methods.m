% evaluate methods to solve the Burgers equation
% Volterra polynomial vs. fractional steps method
%
% author: Martin Schiffner
% date: 2009-04-14
% modified: 2020-05-07

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% physical parameters (distilled water, atmospheric pressure, 20 °C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% independent
%--------------------------------------------------------------------------
mu = 1.002 * 10^-3;        % shear viscosity (Pa s)                (from "Acoustics" by A. D. Pierce, p. 514)
mu_B_over_mu = 2.9;        % bulk viscosity over shear viscosity 3.0 - 2.7 (1) (from "Acoustics" by A. D. Pierce, p. 553)
c_0 = 1482.87;             % small-signal sound speed (m / s) (?)  (from http://www.springerlink.com/content/v04n231880311050/fulltext.pdf)
rho_0 = 998;               % ambient mass density (kg / m^3) (?)   (from http://www.efunda.com/materials/common_matl/show_liquid.cfm?matlname=waterdistilled4c)
B_over_A = 5.0;            % measure of nonlinear effects (1)      (from "Nonlinear Acoustics" by Mark F. Hamilton, p. 34)
beta = 1 + B_over_A / 2;   % coefficient of nonlinearity (1) 
gamma = 1.0079;            % ratio of specific heats               (from "Acoustics" by A. D. Pierce, p. 34)
Pr = 7;                    % Prandtl number mu * c_p / kappa       (from "Acoustics" by A. D. Pierce, p. 514)

%--------------------------------------------------------------------------
% dependent
%--------------------------------------------------------------------------
nu = mu / rho_0;                    % kinematic viscosity (m^2 / s)
delta = nu * (4/3 + mu_B_over_mu + (gamma - 1) / Pr);  % sound diffusivity (m^2 / s) (from "Nonlinear Acoustics" by Mark F. Hamilton)
a = delta / (2 * c_0^3);            % factor a in Burgers' equation
b = beta / (rho_0 * c_0^3);         % factor b in Burgers' equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signal processing parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_s = 1e9;  % sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define input waveform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amplitudes = [ 5e5; 7.5e5 ];%[1; 10; 100; 250; 500; 550; 600; 650; 700; 750; 800; 850; 900; 950; 1000] * 10^3;
N_amplitudes = numel(amplitudes);

%compute Gaussian pulse
f_c = 3.5e6;            % center frequency (Hz)
frac_bw = 1;            % fractional bandwidth
atten_bw = 6;           % attenuation at bandwidth (dB)
atten_td = 100;         % attenuation in time-domain which marks end of Gaussian pulse

%compute dependent parameters
omega_c = 2 * pi * f_c;
sigma = frac_bw * omega_c / (2 * sqrt(atten_bw / (10 * log10(exp(1)))));
N_cutoff = ceil(f_s * sqrt(atten_td / (10 * sigma^2) * log(10)));

tau_axis = (-N_cutoff:N_cutoff) / f_s;
N_samples = 2 * N_cutoff + 1;
pressure_input_temp = gauspuls(tau_axis, f_c, 1, -atten_bw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute reference waveforms using fractional steps method (very small step size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% independent parameters
delta_z_ref = 5e-4;     % propagation distance per step (m)
z_stop = 0.25;          % total propagation distance (m)

% dependent parameters
N_steps_prop = floor( z_stop / delta_z_ref );

% allocate memory for results
pressure_reference = cell( N_amplitudes, N_steps_prop + 1 );
N_steps_nonlin = zeros( N_amplitudes, N_steps_prop + 1 );

% iterate amplitudes
for k = 1:N_amplitudes

    % print status
    fprintf( 'processing amplitude %d of %d (%d)\n', k, N_amplitudes, amplitudes( k ) );

    % define initial pressure waveform
    pressure_reference{ k, 1 } = amplitudes( k ) * pressure_input_temp;

    % compute waveforms for each propagation step
    for l = 1:N_steps_prop

        % print status
        fprintf( '\tstep %d of %d...\n', l, N_steps_prop );

        % call fractional steps method
        %[pressure_reference{k, l + 1}, N_steps_nonlin(k, l + 1)] = burgers_frac_steps(pressure_reference{k, l}, f_s, delta_z_ref, a, b, 0);
        pressure_reference{k, l + 1} = fractional_steps.step( pressure_reference{ k, l }, f_s, delta_z_ref, a, b, 1 );

    end % for l = 1:N_steps_prop

end % for k = 1:N_amplitudes

% save results
str_filename = sprintf( 'pressure_reference_%d.mat', z_stop * 100 );
save( str_filename, 'pressure_reference', 'delta_z_ref', 'f_s', 'amplitudes', 'z_stop' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% evaluate waveforms using Volterra system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% independent parameters
delta_z = 5e-3;                                 % propagation distance per step (m)
orders_N = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ];	% orders of Volterra polynomials (1)

% dependent parameters
N_steps_prop = floor( z_stop / delta_z );	% number of required propagation steps
N_orders = numel( orders_N );

% allocate memory for results
pressure_volterra = cell( N_amplitudes, N_orders, N_steps_prop + 1 );

% iterate amplitudes
for k = 1:N_amplitudes

	% print status
	fprintf( 'processing amplitude %d of %d (%d)\n', k, N_amplitudes, amplitudes( k ) );

	% iterate polynomial orders
	for l = 1:N_orders

        % print status
        fprintf( '\tprocessing polynomial degree %d of %d (%d)\n', l, N_orders, orders_N( l ) );

        % define initial pressure waveform
        pressure_volterra{ k, l, 1 } = amplitudes( k ) * pressure_input_temp;

        % compute waveforms for each propagation step
        for m = 1:N_steps_prop

            % call Volterra method
            pressure_volterra{ k, l, m + 1 } = volterra.polynomial( pressure_volterra{ k, l, m }, f_s, delta_z, a, b, orders_N( l ) );

        end % for m = 1:N_steps_prop

    end % for l = 1:N_orders

end % for k = 1:N_amplitudes

% save results
str_filename = sprintf( 'pressure_volterra_%.1f.mat', delta_z * 1000 );
save( str_filename, 'pressure_volterra', 'orders_N', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% evaluate waveforms using fractional steps method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allocate memory for results
pressure_frac_steps = cell(N_amplitudes, N_steps_prop + 1);

% iterate amplitudes
for k = 1:N_amplitudes 

	% print status
	fprintf( 'processing amplitude %d of %d (%d)\n', k, N_amplitudes, amplitudes( k ) );

	% define initial pressure waveform
	pressure_frac_steps{k,1} = amplitudes(k) * pressure_input_temp;

	% compute waveforms for each propagation step
    for l = 1:N_steps_prop

        % print status
        fprintf( '\tstep %d of %d...\n', l, N_steps_prop );

        % call fractional steps method
%         [ pressure_reference{ k, l + 1 }, N_steps_nonlin( k, l + 1 ) ] = burgers_frac_steps( pressure_reference{ k, l }, f_s, delta_z_ref, a, b, 0 );
        pressure_frac_steps{ k, l + 1 } = fractional_steps.step( pressure_frac_steps{ k, l }, f_s, delta_z, a, b, 1 );

    end % for l = 1:N_steps_prop

end % for k = 1:N_amplitudes 

% save results
str_filename = sprintf( 'pressure_frac_steps_%.1f.mat', delta_z * 1000 );
save( str_filename, 'pressure_frac_steps', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute error metrics of both approaches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error metrics
l2_norm = @(x) sqrt(sum(x.^2));
linfty_norm = @(x) max(abs(x));
correlation_coefficient = @(x,y) xcov(x, y, 0) / (N_samples * sqrt(var(x) * var(y)));
rel_error_l2 = @(p_ref,p_eval) l2_norm(p_ref-p_eval) / l2_norm(p_ref);

% allocate memory for results
error_rel_norm_volterra = zeros( N_amplitudes, N_orders, N_steps_prop );
error_rel_norm_frac_steps = zeros( N_amplitudes, N_steps_prop );
error_rel_max_volterra = zeros( N_amplitudes, N_orders, N_steps_prop );
error_rel_max_frac_steps = zeros( N_amplitudes, N_steps_prop );
coefficient_correlation_volterra = zeros( N_amplitudes, N_orders, N_steps_prop );
coefficient_correlation_frac_steps = zeros( N_amplitudes, N_steps_prop );

% iterate amplitudes
for k = 1:N_amplitudes

    % iterate polynomial orders
    for l = 1:N_orders

        % iterate steps
        for m = 1:N_steps_prop

            pressure_comparison = pressure_reference{ k, m * round( delta_z / delta_z_ref ) + 1 };
            error = pressure_comparison - pressure_volterra{ k, l, m + 1 };
            error_rel_norm_volterra( k, l, m ) = l2_norm( error ) / l2_norm( pressure_comparison );
            error_rel_max_volterra( k, l, m ) = linfty_norm( error ) / linfty_norm( pressure_comparison );
            coefficient_correlation_volterra( k, l, m ) =  correlation_coefficient( pressure_volterra{ k, l, m + 1 }, pressure_comparison );

        end
    end

    % fractional steps method
    for m = 1:N_steps_prop

           pressure_comparison = pressure_reference{ k, m * round( delta_z / delta_z_ref ) + 1 };
           error = pressure_comparison - pressure_frac_steps{ k, m + 1 };
           error_rel_norm_frac_steps( k, m ) = l2_norm( error ) / l2_norm( pressure_comparison );
           error_rel_max_frac_steps( k, m ) = linfty_norm( error ) / linfty_norm( pressure_comparison );
           coefficient_correlation_frac_steps( k, m ) =  correlation_coefficient( pressure_frac_steps{ k, m + 1 }, pressure_comparison );
    end  
end

% save results
str_filename = sprintf( 'errors_%.1f.mat', delta_z * 1000 );
save( str_filename, 'error_rel_norm_volterra', 'error_rel_norm_frac_steps', 'coefficient_correlation_volterra', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% plot error vs. propagation steps for both methods
%--------------------------------------------------------------------------
str_legend = cell(1, N_orders + 1);
for l = 1:N_orders

    str_legend{1, l} = sprintf('M = %d', orders_N(l));
end
str_legend{1, N_orders + 1} = 'frac steps';

for k = 1:N_amplitudes

    %----------------------------------------------------------------------
    % relative error
    %----------------------------------------------------------------------
    figure( 2 * k - 1 );
    hdl_line = semilogy( (1:N_steps_prop), reshape( error_rel_norm_volterra( k, :, : ), N_orders, N_steps_prop )', (1:N_steps_prop), error_rel_norm_frac_steps( k, : ) );
    set( hdl_line( end ), 'LineStyle', '--' );
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [10, 10, 15.5, 9]);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [1.5, 1, 10, 7]);
    set(gca, 'XGrid', 'on');
    set(gca, 'YGrid', 'on');
    
    title( sprintf( 'max pressure: %g kPa, \\Delta z = %d mm', amplitudes( k ) / 1e3, delta_z * 1000 ) );
    xlabel( '# propagation steps (1)' );
    ylabel( 'relative error' );

	hdl_lgnd = legend( str_legend, 'Location', 'eastoutside' );
% 	leg_pos = get( hdl_lgnd, 'Position' );
%     set( hdl_lgnd, 'Position', [ 12, 8 - leg_pos( 4 ), leg_pos( 3 ), leg_pos( 4 ) ] );

    % save figure to file
%     str_filename_emf = sprintf('errors_deltaz_%d_ampl_%d.emf', delta_z * 1000, amplitudes(k));
%     print(gcf, '-r600', '-dmeta', str_filename_emf);
    
    %----------------------------------------------------------------------
    % correlation coefficient
    %----------------------------------------------------------------------
    figure( 2*k );
    hdl_line = semilogy( (1:N_steps_prop), reshape( coefficient_correlation_volterra( k, :, : ), N_orders, N_steps_prop )', (1:N_steps_prop), coefficient_correlation_frac_steps( k, : ) );
    set( hdl_line( end ), 'LineStyle', '--' );

    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [10, 10, 15.5, 9]);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [1.5, 1, 10, 7]);
    set(gca, 'XGrid', 'on');
    set(gca, 'YGrid', 'on');
    
    title(sprintf('max pressure: %g kPa, \\Delta z = %d mm', amplitudes( k ) / 1e3, delta_z * 1000));
    xlabel('# propagation steps (1)');
    ylabel('correlation coefficient');
    
    hdl_lgnd = legend( str_legend, 'Location', 'eastoutside' );
%     leg_pos = get(hdl_lgnd, 'Position');
%     set(hdl_lgnd, 'Position', [12, 8 - leg_pos(4), leg_pos(3), leg_pos(4)]);
    
end

%--------------------------------------------------------------------------
% plot error vs. order of Volterra polynomial
%--------------------------------------------------------------------------
str_legend_2 = cell(1, N_amplitudes);
for l = 1:N_amplitudes

    str_legend_2{1, l} = sprintf('%d kPa', amplitudes(l) / 1e3);
end

figure( 2*N_amplitudes + 1 );
hdl_line = semilogy( orders_N, error_rel_norm_volterra( :, :, N_steps_prop )', orders_N, repmat( error_rel_norm_frac_steps( :, N_steps_prop ), [ 1, numel( orders_N ) ] ) );
set( hdl_line( (end - N_amplitudes + 1):end ), 'LineStyle', '--' );

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 15.5, 9]);
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1.5, 1, 10, 7]);
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');

title(sprintf('z = %d mm', N_steps_prop * delta_z * 1000));
xlabel( 'Degree of Volterra polynomial (1)' );
ylabel( 'Relative RMSE (1)' );
    
hdl_lgnd = legend( str_legend_2, 'Location', 'eastoutside' );
% leg_pos = get(hdl_lgnd, 'Position');
% set(hdl_lgnd, 'Position', [12, 8 - leg_pos(4), leg_pos(3), leg_pos(4)]);

%save figure to file
% str_filename_emf = sprintf('errors_deltaz_%d_final.emf', delta_z * 1000);
% print(gcf, '-r600', '-dmeta', str_filename_emf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create movie (high pressure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index_figure = 2 * N_amplitudes + 2;
index_amplitude = 2;    % select amplitude
index_order = 10;       % order of Volterra polynomial

indices_steps_draw = (1:1:N_steps_prop + 1);
N_steps_draw = numel( indices_steps_draw );

N_dft = 8192;
index_shift = ceil( N_dft / 2 );
f_axis = f_s * (0:(index_shift - 1)) / N_dft;

%--------------------------------------------------------------------------
% geometrical parameters
%--------------------------------------------------------------------------
scale = 1;
N_plots_x = 2;
N_plots_y = 1;
width_axis_x = 8 * scale;
width_axis_y = 7 * scale;
spacing_axis_x = 1.5 * ones( 1, N_plots_x - 1 );
spacing_axis_y = 0.75;
axis_x_offset = 1.5;
axis_y_offset = 1.5;
width_cb_x = 0.4;
delta_cb = 0.3;
delta_lgd = 0.4;
height_lgd = 0.75;

text_size = 10;

size_legend = size_label;

%--------------------------------------------------------------------------
% axis limits and ticks
%--------------------------------------------------------------------------
limits_tau = [ -0.6, 0.6 ];
limits_p = [ -0.5, 1.1 ];

limits_f = [ 0, 20 ];
limits_dB = [ -40, 0 ];

tau_tick = (-0.6:0.1:0.6);
p_tick = (-0.5:0.25:1.1);

f_tick = (0:1:20);
dB_tick = (-40:5:0);

ticks_font = 'Times';
ticks_size = text_size;
ticks_color = 'k';

%--------------------------------------------------------------------------
% axis labels
%--------------------------------------------------------------------------
label_font = 'Times';
label_size = text_size + 2;
label_color = 'k';

label_tau_str = 'Retarded time (us)';
label_p_str = 'Acoustic pressure (MPa)';

label_f_str = 'Normalized frequency (1)';
label_dB_str = 'Spectrum (dB)';

%--------------------------------------------------------------------------
% titles
%--------------------------------------------------------------------------
title_font = 'Times';
title_size = text_size + 4;
title_color = 'k';

title_factor_shift_x = 0;
title_factor_shift_y = -0.075;

%--------------------------------------------------------------------------
% dependent parameters
%--------------------------------------------------------------------------
pos_x = zeros( 1, N_plots_x );
pos_y = zeros( 1, N_plots_y );

for index_x = 1:N_plots_x
    pos_x( index_x ) = axis_x_offset + ( index_x - 1 ) * width_axis_x + sum( spacing_axis_x( 1:( index_x - 1 ) ) );
end
for index_y = 1:N_plots_y
    pos_y( index_y ) = axis_y_offset + ( N_plots_y - index_y ) * ( width_axis_y + spacing_axis_y );
end

width_cb_y = pos_y(1) + width_axis_y - pos_y(end);

pos_x_cb = pos_x(end) + width_axis_x + delta_cb;

paper_size_x = pos_x( end ) + width_axis_x + delta_cb + width_cb_x + 0.5;
paper_size_y = pos_y( 1 ) + width_axis_y + height_lgd + 0.5;

% create a video writer object for the output video file and open the object for writing
str_filename = 'burgers_propagation_hp.gif';
% video = VideoWriter( 'burgers_propagation_hp.avi', 'Motion JPEG AVI' );
% video.FrameRate = 10;
% open( video );

%--------------------------------------------------------------------------
% create figure
%--------------------------------------------------------------------------
figure( index_figure );

set( gcf, 'Units', 'centimeters');
set( gcf, 'Position', [ 10, 10, paper_size_x, paper_size_y ] );
set( gcf, 'PaperSize', [ paper_size_x, paper_size_y ] );
set( gcf, 'PaperPositionMode', 'auto');
set(gcf, 'NextPlot', 'replacechildren');

axes_hdl = zeros( 1, N_plots_x * N_plots_y );

for k = 1:N_steps_draw

	%----------------------------------------------------------------------
	% plot waveform
	%----------------------------------------------------------------------
	axes_hdl( 1 ) = subplot( 1, 2, 1 );
    plot( tau_axis * 1e6, pressure_volterra{ index_amplitude, 1, indices_steps_draw( k ) } / 1e6, ...
          tau_axis * 1e6, pressure_volterra{ index_amplitude, index_order, indices_steps_draw( k ) } / 1e6 );

	%----------------------------------------------------------------------
	% plot spectrum
	%----------------------------------------------------------------------
	axes_hdl( 2 ) = subplot( 1, 2, 2 );
	pressure_reference_dft = abs( fft( pressure_volterra{ index_amplitude, 1, indices_steps_draw( k ) }, N_dft ) );
	pressure_reference_dft_dB = 20 * log10( pressure_reference_dft / max( pressure_reference_dft ) );
	pressure_volterra_dft = abs( fft( pressure_volterra{ index_amplitude, index_order, indices_steps_draw( k ) }, N_dft ) );
	pressure_volterra_dft_dB = 20 * log10( pressure_volterra_dft / max(pressure_volterra_dft ) );

    plot( f_axis / f_c, pressure_reference_dft_dB( 1:index_shift ), f_axis / f_c, pressure_volterra_dft_dB( 1:index_shift ) );

    %----------------------------------------------------------------------
    % legend
    %----------------------------------------------------------------------
	leg_hdl = legend( { 'linear propagation', 'nonlinear propagation' }, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Units', 'centimeters' );
    set( leg_hdl, 'FontSize', size_legend );
	set( leg_hdl, 'Position', [ pos_x( 2 ), pos_y( 1 ) + width_axis_y + delta_lgd, width_axis_x, height_lgd ] );

    %----------------------------------------------------------------------
    % format plots
    %----------------------------------------------------------------------
    for index_axis = 1:numel( axes_hdl )

        % size and position
        set( axes_hdl( index_axis ), 'Units', 'centimeters' );
        set( axes_hdl( index_axis ), 'FontName', ticks_font, 'FontSize', ticks_size );

        if index_axis == 1
            set( axes_hdl( index_axis ), 'Position', [ pos_x( index_axis ), pos_y( 1 ), width_axis_x, width_axis_y ] );
        else
            set( axes_hdl( index_axis ), 'Position', [ pos_x( index_axis ), pos_y( 1 ), width_axis_x, width_axis_y ] );
        end

        % axes limits and ticks
        if index_axis == 1
            xlim( axes_hdl( index_axis ), limits_tau );
            ylim( axes_hdl( index_axis ), limits_p );
            set( axes_hdl( index_axis ), 'XTick', tau_tick );
            set( axes_hdl( index_axis ), 'YTick', p_tick );
        else
            xlim( axes_hdl( index_axis ), limits_f );
            ylim( axes_hdl( index_axis ), limits_dB );
            set( axes_hdl( index_axis ), 'XTick', f_tick );
            set( axes_hdl( index_axis ), 'YTick', dB_tick );
        end
        set( axes_hdl( index_axis ), 'XGrid', 'on' );
        set( axes_hdl( index_axis ), 'YGrid', 'on' );

        % labels
        if index_axis == 1
            xlbl_hdl = xlabel( axes_hdl( index_axis ), label_tau_str, 'FontName', label_font, 'FontSize', label_size, 'Color', label_color );
            xlbl_pos = get( xlbl_hdl, 'Position' );
            set( xlbl_hdl, 'Position', xlbl_pos + [0, 0, 0] );
        else
            xlbl_hdl = xlabel( axes_hdl( index_axis ), label_f_str, 'FontName', label_font, 'FontSize', label_size, 'Color', label_color );
%             set(xlbl_hdl, 'Position', xlbl_pos + [0, 0, 0]);
        end

        if index_axis == 1
            ylbl_hdl = ylabel( axes_hdl( index_axis ), label_p_str, 'FontName', label_font, 'FontSize', label_size, 'Color', label_color );
            ylbl_pos = get( ylbl_hdl, 'Position' );
            set(ylbl_hdl, 'Position', ylbl_pos + [0, 0, 0]);
        else
            ylbl_hdl = ylabel( axes_hdl( index_axis ), label_dB_str, 'FontName', label_font, 'FontSize', label_size, 'Color', label_color );
%             set(ylbl_hdl, 'Position', ylbl_pos + [0, 0, 0]);
        end

        % titles
        if index_axis == 1
            set( gcf, 'currentaxes', axes_hdl( index_axis ) );
            tit_pos_x = limits_tau( 1 ) + title_factor_shift_x * diff( limits_tau );
            tit_pos_y = limits_p( 2 ) - title_factor_shift_y * diff( limits_p );
            tit_hdl = text( tit_pos_x, tit_pos_y, sprintf( 'Propagation distance: %.1f cm', ( indices_steps_draw( k ) - 1 ) * delta_z * 1e2 ), 'Units', 'data', 'FontName', title_font, 'FontSize', title_size, 'Interpreter', 'none', 'Color', title_color, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left' );
        end

    end % for index_axis = 1:numel( axes_hdl )

    %----------------------------------------------------------------------
    % write figure to movie
    %----------------------------------------------------------------------
	% write current figure as frame
    frame = getframe( gcf );
% 	writeVideo( video, frame );
    image = frame2im( frame );
	[ imind, cm ] = rgb2ind( image, 256 );

	% write to the GIF File
	if k == 1
        imwrite( imind, cm, str_filename, 'gif', 'Loopcount', inf);
    else
        imwrite( imind, cm, str_filename, 'gif', 'WriteMode', 'append');
    end

end % for k = 1:N_steps_draw

% close video object
% close( video );
