%evaluate methods to solve the Burgers equation
%Volterra polynomial vs. fractional steps solution
%author: Martin Schiffner
%date: 2009-04-14

%% clear workspace
clear;
close all;
clc;

%% physical parameters (distilled water, atmospheric pressure, 20 �C)
mu = 1.002 * 10^-3;        %shear viscosity (Pa s)                (from "Acoustics" by A. D. Pierce, p. 514)
mu_B_over_mu = 2.9;        %bulk viscosity over shear viscosity 3.0 - 2.7 (1) (from "Acoustics" by A. D. Pierce, p. 553)
c_0 = 1482.87;             %small-signal sound speed (m / s) (?)  (from http://www.springerlink.com/content/v04n231880311050/fulltext.pdf)
rho_0 = 998;               %ambient mass density (kg / m^3) (?)   (from http://www.efunda.com/materials/common_matl/show_liquid.cfm?matlname=waterdistilled4c)
B_over_A = 5.0;            %measure of nonlinear effects (1)      (from "Nonlinear Acoustics" by Mark F. Hamilton, p. 34)
beta = 1 + B_over_A / 2;   %coefficient of nonlinearity (1) 
gamma = 1.0079;            %ratio of specific heats               (from "Acoustics" by A. D. Pierce, p. 34)
Pr = 7;                    %Prandtl number mu * c_p / kappa       (from "Acoustics" by A. D. Pierce, p. 514)

%% dependent physical parameters
nu = mu / rho_0;                    %kinematic viscosity (m^2 / s)
delta = nu * (4/3 + mu_B_over_mu + (gamma - 1) / Pr);  %sound diffusivity (m^2 / s) (from "Nonlinear Acoustics" by Mark F. Hamilton)
a = delta / (2 * c_0^3);            %factor a in Burgers' equation
b = beta / (rho_0 * c_0^3);         %factor b in Burgers' equation

%% parameters for signal processing
f_s = 10^9;

%% define input waveform
amplitudes = [1; 10; 100; 250; 500; 550; 600; 650; 700; 750; 800; 850; 900; 950; 1000] * 10^3;
N_amplitudes = numel(amplitudes);

%compute Gaussian pulse
f_c = 3.5 * 10^6;       % center frequency (Hz)
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

%% compute reference waveforms using fractional steps approach

delta_z_ref = 0.0005;   % propagation distance per step (m)
z_stop = 0.25;          % total propagation distance (m)
N_steps_prop = z_stop / delta_z_ref;

%allocate memory for results
pressure_reference = cell(N_amplitudes, N_steps_prop + 1);
N_steps_nonlin = zeros(N_amplitudes, N_steps_prop + 1);

for k = 1:N_amplitudes 
    
    %print status
    disp(sprintf('processing amplitude %d of %d (%d)', k, N_amplitudes, amplitudes(k)));
    
    %define initial pressure waveform
    pressure_reference{k,1} = amplitudes(k) * pressure_input_temp;
    
    %compute waveforms for each propagation step
    for l = 1:N_steps_prop
        
        disp(sprintf('\tstep %d of %d...', l, N_steps_prop));
        %[pressure_reference{k, l + 1}, N_steps_nonlin(k, l + 1)] = burgers_frac_steps(pressure_reference{k, l}, f_s, delta_z_ref, a, b, 0);
        pressure_reference{k, l + 1} = burgers_frac_steps_m(pressure_reference{k, l}, f_s, delta_z_ref, a, b, 1);
    end
end

%save results
save(sprintf('pressure_reference_%d.mat', z_stop * 100), 'pressure_reference', 'delta_z_ref', 'f_s', 'amplitudes', 'z_stop');

%% evaluate waveforms using Volterra system

delta_z = 0.005;    %propagation distance per step (m)
N_steps_prop = floor(z_stop / delta_z);     %number of required propagation steps
orders_N = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
N_orders = numel(orders_N);

%allocate memory for results
pressure_volterra = cell(N_amplitudes, N_orders, N_steps_prop + 1);

for k = 1:N_amplitudes
    
    %print status
    disp(sprintf('processing amplitude %d of %d (%d),', k, N_amplitudes, amplitudes(k)));
    
    for l = 1:N_orders
        
        %print status
        disp(sprintf('\tprocessing polynomial degree %d of %d (%d),', l, N_orders, orders_N(l)));
        
        %define initial pressure waveform
        pressure_volterra{k, l, 1} = amplitudes(k) * pressure_input_temp;
        
        %compute waveforms for each propagation step
        for m = 1:N_steps_prop
            
            %disp(sprintf('\t\tstep %d of %d...', m, N_steps_prop));
            pressure_volterra{k, l, m + 1} = burgers_volterra(pressure_volterra{k, l, m}, f_s, delta_z, a, b, orders_N(l));
        end
    end
end

%save results
save(sprintf('pressure_volterra_%.1f.mat', delta_z * 1000), 'pressure_volterra', 'orders_N', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop');

%% evaluate waveforms using fractional steps approach

%allocate memory for results
pressure_frac_steps = cell(N_amplitudes, N_steps_prop + 1);

for k = 1:N_amplitudes 
    
    %print status
    disp(sprintf('processing amplitude %d of %d (%d)', k, N_amplitudes, amplitudes(k)));
    
    %define initial pressure waveform
    pressure_frac_steps{k,1} = amplitudes(k) * pressure_input_temp;
    
    %compute waveforms for each propagation step
    for l = 1:N_steps_prop
    
        disp(sprintf('\tstep %d of %d...', l, N_steps_prop));
        %[pressure_reference{k, l + 1}, N_steps_nonlin(k, l + 1)] = burgers_frac_steps(pressure_reference{k, l}, f_s, delta_z_ref, a, b, 0);
        pressure_frac_steps{k, l + 1} = burgers_frac_steps_m(pressure_frac_steps{k, l}, f_s, delta_z, a, b, 1);
    end
end

%save results
save(sprintf('pressure_frac_steps_%.1f.mat', delta_z * 1000), 'pressure_frac_steps', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop');

%% compute error metrics of both approaches

l2_norm = @(x) sqrt(sum(x.^2));
linfty_norm = @(x) max(abs(x));
correlation_coefficient = @(x,y) xcov(x, y, 0) / (N_samples * sqrt(var(x) * var(y)));
rel_error_l2 = @(p_ref,p_eval) l2_norm(p_ref-p_eval) / l2_norm(p_ref);

%allocate memory for results
error_rel_norm_volterra = zeros(N_amplitudes, N_orders, N_steps_prop);
error_rel_norm_frac_steps = zeros(N_amplitudes, N_steps_prop);
error_rel_max_volterra = zeros(N_amplitudes, N_orders, N_steps_prop);
error_rel_max_frac_steps = zeros(N_amplitudes, N_steps_prop);
coefficient_correlation_volterra = zeros(N_amplitudes, N_orders, N_steps_prop);

for k = 1:N_amplitudes 
    
    %error caused by Volterra polynomial
    for l = 1:N_orders
        
        for m = 1:N_steps_prop
    
            pressure_comparison = pressure_reference{k, m * floor(delta_z / delta_z_ref) + 1};
            error = pressure_comparison - pressure_volterra{k, l, m + 1};
            error_rel_norm_volterra(k, l, m) = l2_norm(error) / l2_norm(pressure_comparison);
            error_rel_max_volterra(k, l, m) = linfty_norm(error) / linfty_norm(pressure_comparison);
            coefficient_correlation_volterra(k, l, m) =  correlation_coefficient(pressure_volterra{k, l, m + 1}, pressure_comparison);
        end
    end
    
    %error caused by fractional steps approach
    for m = 1:N_steps_prop
     
           pressure_comparison = pressure_reference{k, m * floor(delta_z / delta_z_ref) + 1};
           error = pressure_comparison - pressure_frac_steps{k, m + 1};
           error_rel_norm_frac_steps(k, m) = l2_norm(error) / l2_norm(pressure_comparison);
           error_rel_max_frac_steps(k, m) = linfty_norm(error) / linfty_norm(pressure_comparison);
    end  
end

%save results
save(sprintf('errors_%.1f.mat', delta_z * 1000), 'error_rel_norm_volterra', 'error_rel_norm_frac_steps', 'coefficient_correlation_volterra', 'delta_z', 'f_s', 'amplitudes', 'z_stop', 'N_steps_prop');

%% plot results

%plot error vs. propagation steps for both methods
str_legend = cell(1, N_orders + 1);
for l = 1:N_orders

    str_legend{1, l} = sprintf('M = %d', orders_N(l));
end
str_legend{1, N_orders + 1} = 'frac steps';

for k = 1:N_amplitudes
    
    figure(2*k-1);
    %hdl_line = semilogy((1:N_steps_prop), reshape(error_rel_norm_volterra(k,:,:), N_orders, N_steps_prop)', (1:N_steps_prop), error_rel_norm_frac_steps(k,:));   
    hdl_line = semilogy((1:N_steps_prop), reshape(error_rel_norm_volterra(k,:,:), N_orders, N_steps_prop)');   
    %set(hdl_line(1:end-1), 'Marker', 's', 'MarkerSize', 5);
    set(hdl_line(end), 'LineStyle', '--');
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [10, 10, 15.5, 9]);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [1.5, 1, 10, 7]);
    
    title(sprintf('max pressure: %g Pa, \\Delta z = %d mm', amplitudes(k), delta_z * 1000));
    xlabel('# propagation steps');
    ylabel('relative error');
    
    hdl_lgnd = legend(str_legend);
    leg_pos = get(hdl_lgnd, 'Position');
    set(hdl_lgnd, 'Position', [12, 8 - leg_pos(4), leg_pos(3), leg_pos(4)]);
    
    %save figure to file
    str_filename_emf = sprintf('errors_deltaz_%d_ampl_%d.emf', delta_z * 1000, amplitudes(k));
    print(gcf, '-r600', '-dmeta', str_filename_emf);
    
    figure(2*k);
    hdl_line = semilogy((1:N_steps_prop), reshape(coefficient_correlation_volterra(k,:,:), N_orders, N_steps_prop)');
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [10, 10, 15.5, 9]);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [1.5, 1, 10, 7]);
    
    title(sprintf('max pressure: %g Pa, \\Delta z = %d mm', amplitudes(k), delta_z * 1000));
    xlabel('# propagation steps');
    ylabel('correlation coefficient');
    
    hdl_lgnd = legend(str_legend);
    leg_pos = get(hdl_lgnd, 'Position');
    set(hdl_lgnd, 'Position', [12, 8 - leg_pos(4), leg_pos(3), leg_pos(4)]);
end

%plot error vs. order of Volterra polynomial

str_legend_2 = cell(1, N_amplitudes);
for l = 1:N_amplitudes

    str_legend_2{1, l} = sprintf('%d Pa', amplitudes(l));
end

figure(N_amplitudes + 1);
hdl_line = semilogy(orders_N, error_rel_norm_volterra(:,:,N_steps_prop)', orders_N(end) + 1, error_rel_norm_frac_steps(:,N_steps_prop));
set(hdl_line((end - N_amplitudes + 1):end), 'Marker', 's', 'MarkerSize', 5);

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 15.5, 9]);
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1.5, 1, 10, 7]);

title(sprintf('z = %d mm', N_steps_prop * delta_z * 1000));
xlabel('degree of Volterra polynomial');
ylabel('relative error');
    
hdl_lgnd = legend(str_legend_2);
leg_pos = get(hdl_lgnd, 'Position');
set(hdl_lgnd, 'Position', [12, 8 - leg_pos(4), leg_pos(3), leg_pos(4)]);

%save figure to file
str_filename_emf = sprintf('errors_deltaz_%d_final.emf', delta_z * 1000);
print(gcf, '-r600', '-dmeta', str_filename_emf);

%% create movie (high pressure)
draw_steps = (1:2:N_steps_prop + 1);
N_draw_steps = numel(draw_steps);

index_nyquist = ceil(numel(tau_axis) / 2);
f_axis = f_s * (0:(index_nyquist - 1)) / numel(tau_axis);

F = moviein(N_draw_steps);
figure(1);
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 20, 13]);
set(gcf, 'PaperSize', [10, 10]);
set(gcf, 'NextPlot', 'replacechildren');

for k = 1:N_draw_steps
    
    %plot waveform
    subplot(1,2,1);   
    plot(gca, tau_axis * 10^6, pressure_reference{7,draw_steps(k)} / 10^6, tau_axis * 10^6, pressure_volterra{7,1,draw_steps(k)} / 10^6);
    
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [1.5, 1.5, 8, 11]);
	set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    
    xlim([-0.6,0.6]);
    set(gca, 'XTick', (-0.6:0.1:0.6));
    set(gca, 'XGrid', 'on');
    xlbl_hdl = xlabel('$\tau$ / $\mu$s', 'Interpreter', 'LaTeX', 'FontSize', 16);
    xlbl_pos = get(xlbl_hdl, 'Position');
    set(xlbl_hdl, 'Position', xlbl_pos + [0, 0, 0]);
    ylim([-0.5, 1.1]);
    set(gca, 'YGrid', 'on');
    ylbl_hdl = ylabel('p / MPa', 'Interpreter', 'LaTeX', 'FontSize', 16);
    ylbl_pos = get(ylbl_hdl, 'Position');
    set(ylbl_hdl, 'Position', ylbl_pos + [0.02, 0, 0]);
   
    text(-0.2, 1.05, sprintf('$z = %.1f$ cm', (draw_steps(k) - 1) * delta_z * 100), 'Interpreter', 'LaTeX', 'FontSize', 16);
    
    %plot spectrum
    subplot(1,2,2);
    pressure_reference_dft = abs(fft(pressure_reference{7,draw_steps(k)}));
    pressure_reference_dft_dB = 20 * log10(pressure_reference_dft / max(pressure_reference_dft));
    pressure_volterra_dft = abs(fft(pressure_volterra{7,1,draw_steps(k)})) / 10^6;
    pressure_volterra_dft_dB = 20 * log10(pressure_volterra_dft / max(pressure_volterra_dft));
    
    plot(gca, f_axis / 10^6, pressure_reference_dft_dB(1:index_nyquist), f_axis / 10^6, pressure_volterra_dft_dB(1:index_nyquist));
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [11, 1.5, 8, 11]);
	set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    
    
    xlim([0, 30]);
    set(gca, 'XGrid', 'on');
    xlbl_hdl = xlabel('$f$ / MHz', 'Interpreter', 'LaTeX', 'FontSize', 16);
    xlbl_pos = get(xlbl_hdl, 'Position');
    set(xlbl_hdl, 'Position', xlbl_pos + [0, 0, 0]);
    ylim([0, 30]);
    
    leg_hdl = legend('nichtlinear', 'linear');
    set(leg_hdl, 'Interpreter', 'LaTeX', 'FontSize', 16);
    leg_pos = get(leg_hdl, 'Position');
    set(leg_hdl, 'Position', [9.5 - leg_pos(3), 12.5 - leg_pos(4), leg_pos(3), leg_pos(4)]);
    
    F(k) = getframe(gcf);
end

movie2avi(F, 'C:\temp\burgers_comparison_hp.avi');

%% create movie (low pressure)
draw_steps = (1:2:N_steps_prop + 1);
N_draw_steps = numel(draw_steps);

F = moviein(N_draw_steps);
figure(1);
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 18, 13]);
set(gcf, 'PaperSize', [10, 10]);
set(gcf, 'NextPlot', 'replacechildren');

set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1.5, 1.5, 16, 11]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

for k = 1:N_draw_steps
    
    %plot waveform   
    plot(gca, tau_axis * 10^6, pressure_reference{3,draw_steps(k)} / 10^3, tau_axis * 10^6, pressure_volterra{3,1,draw_steps(k)} / 10^3);
    xlim([-0.6,0.6]);
    set(gca, 'XTick', (-0.6:0.1:0.6));
    set(gca, 'XGrid', 'on');
    xlbl_hdl = xlabel('$\tau$ / $\mu$s', 'Interpreter', 'LaTeX', 'FontSize', 16);
    xlbl_pos = get(xlbl_hdl, 'Position');
    set(xlbl_hdl, 'Position', xlbl_pos + [0, 0, 0]);
    ylim([-50, 110]);
    set(gca, 'YGrid', 'on');
    ylbl_hdl = ylabel('p / kPa', 'Interpreter', 'LaTeX', 'FontSize', 16);
    ylbl_pos = get(ylbl_hdl, 'Position');
    set(ylbl_hdl, 'Position', ylbl_pos + [0.02, 0, 0]);
        
    leg_hdl = legend('nichtlinear', 'linear');
    set(leg_hdl, 'Interpreter', 'LaTeX', 'FontSize', 16);
    leg_pos = get(leg_hdl, 'Position');
    set(leg_hdl, 'Position', [17.5 - leg_pos(3), 12.5 - leg_pos(4), leg_pos(3), leg_pos(4)]);

    text(0.36, 75, sprintf('$z = %.1f$ cm', (draw_steps(k) - 1) * delta_z * 100), 'Interpreter', 'LaTeX', 'FontSize', 16);  
   
    F(k) = getframe(gcf);
end

movie2avi(F, 'C:\temp\burgers_comparison_lp.avi');

%% create figures for IUS paper
index_amplitude = 5;      %amplitude used to plot waveforms
max_distance = z_stop;    %maximum propagation distance (m)
N_steps_max_ref = max_distance / delta_z_ref + 1;
delta_z_eval = [2.5, 5, 7.5, 10, 12.5];
%numerische Probleme aufgrund vieler Ausbreitungsschritte bei 1mm

%initial waveform in time-domain (normalized)
figure(1)
plot(tau_axis * 10^9, pressure_input_temp);

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 4.3, 3.5]);
set(gcf, 'PaperSize', [4.3, 3.5]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1, 0.8, 3, 2.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 6);

xlim ([-400, 400]);
set(gca, 'XTick', (-400:200:400));
set(gca, 'XMinorTick', 'off');
set(gca, 'XGrid', 'on');

ylim ([-0.5, 1]);
set (gca, 'YTick', (-0.5:0.25:1));
set (gca, 'YMinorTick', 'off');
set (gca, 'YGrid', 'on');
set(gca, 'GridLineStyle', '-.');

xl_hdl = xlabel ('$\tau$ / ns', 'Interpreter', 'LaTeX', 'FontSize', 8);
xl_pos = get (xl_hdl, 'Position');
set (xl_hdl, 'Position', xl_pos + [0, 0.1, 0]);

ylbl_hdl = ylabel ('$p_{0}(\tau) / \left\|p_{0}\right\|_{\infty}$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
ylbl_pos = get(ylbl_hdl, 'Position');
set (ylbl_hdl, 'Position', ylbl_pos + [0.55, 0, 0]);

print('-r600', '-dpdf', 'pressure_initial.pdf');

%waveforms at maximum propagation distance
figure(2);
line_hdl = plot(tau_axis * 10^9, pressure_reference{index_amplitude, N_steps_max_ref} / amplitudes(index_amplitude),...
                tau_axis * 10^9, pressure_volterra{index_amplitude, 1, N_steps_max} / amplitudes(index_amplitude),...
                tau_axis * 10^9, pressure_volterra{index_amplitude, 2, N_steps_max} / amplitudes(index_amplitude),...
                tau_axis * 10^9, pressure_volterra{index_amplitude, 3, N_steps_max} / amplitudes(index_amplitude));
set(line_hdl(1), 'LineWidth', 0.8, 'LineStyle', '--', 'Color', 'r');
set(line_hdl(2), 'LineWidth', 0.4, 'LineStyle', '-', 'Color', 'c');
set(line_hdl(3), 'LineWidth', 0.4, 'LineStyle', '-','Color', 'g');
set(line_hdl(4), 'LineWidth', 0.4, 'LineStyle', '-', 'Color', 'b');

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 4.3, 3.5]);
set(gcf, 'PaperSize', [4.3, 3.5]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1, 0.8, 3, 2.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 6);

xlim ([-400, 400]);
set(gca, 'XTick', (-400:200:400));
set(gca, 'XMinorTick', 'off');
set(gca, 'XGrid', 'on');

ylim([-0.5, 1]);
set(gca, 'YTick', (-0.5:0.25:1));
set(gca, 'YMinorTick', 'off');
set(gca, 'YGrid', 'on');
set(gca, 'GridLineStyle', '-.');

xl_hdl = xlabel ('$\tau$ / ns', 'Interpreter', 'LaTeX', 'FontSize', 8);
xl_pos = get (xl_hdl, 'Position');
set(xl_hdl, 'Position', xl_pos + [0, 0.1, 0]);

ylbl_hdl = ylabel ('$p_{z_{\mathrm{max}}}^{(M_{i})}(\tau) / \left\|p_{0}\right\|_{\infty}$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
ylbl_pos = get(ylbl_hdl, 'Position');
set (ylbl_hdl, 'Position', ylbl_pos + [0.55, 0, 0]);

print('-r600', '-dpdf', 'pressure_propagated.pdf');

%rel. error vs. propagation distance
figure(3);
data = reshape(error_rel_norm_volterra(index_amplitude,:,:), N_orders, N_steps_prop);
line_hdl = semilogy((1:N_steps_prop), data(:, 1:N_steps_prop));
%line_hdl = semilogy((1:N_steps_prop), data(:, 1:N_steps_prop), (1:N_steps_prop), error_rel_norm_frac_steps(index_amplitude,:));

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 8.8, 4.61]);
set(gcf, 'PaperSize', [8.8, 4.61]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1, 0.8, 7.6, 3.8]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 6);

xlim ([0, N_steps_prop]);
set(gca, 'XTick', (0:5:N_steps_prop));
set(gca, 'XMinorTick', 'off');
set(gca, 'XGrid', 'off');

ylim([4*10^-4, 1]);
set(gca, 'YTick', [10^-3, 10^-2, 10^-1]);
%set(gca, 'YTickLabel', ytick_hdl);
set(gca, 'YMinorTick', 'off');
set(gca, 'YGrid', 'on');
set(gca, 'GridLineStyle', '-.');

xl_hdl = xlabel ('$z / \Delta z$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
xl_pos = get (xl_hdl, 'Position');
set(xl_hdl, 'Position', xl_pos + [0, 0.15, 0]);

ylbl_hdl = ylabel ('$\varepsilon^{(M)}(z)$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
ylbl_pos = get(ylbl_hdl, 'Position');
set (ylbl_hdl, 'Position', ylbl_pos + [0.4, 0, 0]);

text(2, 0.5, '$M = 1$', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', 'b');
line([4, 5], [0.38, 0.133], 'Color', 'b');
text(5, 0.04, '$M = 2$', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', [0,0.5,0]);
line([7, 8], [0.03, 0.01347], 'Color', [0,0.5,0]);
text(15, 0.09, '$M = 3$', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', 'r');
line([17, 19], [0.07, 0.01709], 'Color', 'r');
text(21, 0.007, '$M = 4-~10$', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', 'r');
line([21, 20], [0.008, 0.016], 'Color', 'r');

print('-r600', '-dpdf', 'error_rel_l2_vs_z.pdf');
 
%rel. error at max distance vs. degree M (for delta_z = 1, 2, 3, 4, 5, 7.5, 10, 12.5, 15 mm)
%load relative errors
error_rel_norm_volterra_eval = zeros(numel(delta_z_eval), N_orders);
error_rel_norm_frac_steps_eval = zeros(numel(delta_z_eval), 1);

for k = 1:numel(delta_z_eval)

    str_filename = sprintf('errors_%.1f.mat', delta_z_eval(k));
    temp = load(str_filename);
    error_rel_norm_volterra_eval(k,:) = temp.error_rel_norm_volterra(index_amplitude,:,temp.N_steps_prop);
    error_rel_norm_frac_steps_eval(k) = temp.error_rel_norm_frac_steps(index_amplitude,temp.N_steps_prop);
end
figure(4);
line_hdl = semilogy((1:N_orders), error_rel_norm_volterra_eval);
set(line_hdl(:), 'LineWidth', 0.4, 'Marker', '+', 'MarkerSize', 6);
    
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 8.8, 4.43]);
set(gcf, 'PaperSize', [8.8, 4.43]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1, 0.55, 7.6, 3.8]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 6);
    
xlim ([1, 4]);
set(gca, 'XTick', (1:N_orders));
set(gca, 'XMinorTick', 'off');
set(gca, 'XGrid', 'off');

ylim([10^-1, 1]);
set(gca, 'YTick', (0.1:0.1:1));
%set(gca, 'YTickLabel', [0.1, 0.25, 0.5, 0.75, 1]);
set(gca, 'YMinorTick', 'off');
set(gca, 'YGrid', 'on');
set(gca, 'GridLineStyle', '-.');

xl_hdl = xlabel ('$M$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
xl_pos = get (xl_hdl, 'Position');
set(xl_hdl, 'Position', xl_pos + [0, 0.3, 0]);

ylbl_hdl = ylabel ('$\varepsilon^{(M)}(z_{\mathrm{max}})$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
ylbl_pos = get(ylbl_hdl, 'Position');
set (ylbl_hdl, 'Position', ylbl_pos + [0.35, 0, 0]);

text(1.1, 0.12, '$\Delta z = 2.5$ mm', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', 'b');
line([1.4, 1.78], [0.13, 0.2], 'Color', 'b');
text(1.1, 0.22, '$\Delta z = 5$ mm', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', [0,0.5,0]);
line([1.35, 1.6], [0.24, 0.3], 'Color', [0,0.5,0]);
text(1.9, 0.25, '$\Delta z = 7.5$ mm', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', 'r');
line([2.25, 2.15], [0.146, 0.22], 'Color', 'r');
text(1.6, 0.45, '$\Delta z = 10$ mm', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', [0,0.75,0.75]);
line([1.9, 1.68], [0.4, 0.3], 'Color', [0,0.75,0.75]);
text(2.5, 0.5, '$\Delta z = 12.5$ mm', 'Interpreter', 'LaTeX', 'FontSize', 8, 'Color', [0.75,0,0.75]);
line([2.8, 2.5], [0.45, 0.3], 'Color', [0.75,0,0.75]);

print('-r600', '-dpdf', 'error_rel_l2_vs_M.pdf');

close(gcf);

%% create figures for IUS presentation
index_amplitude = 5;      %amplitude used to plot waveforms
max_distance = z_stop;    %maximum propagation distance (m)
N_steps_max_ref = max_distance / delta_z_ref + 1;
delta_z_eval = [2.5, 5, 7.5, 10, 12.5, 15];
%numerische Probleme aufgrund vieler Ausbreitungsschritte bei 1mm

%initial waveform in time-domain (normalized)
figure(1)
plot(tau_axis * 10^9, pressure_input_temp, 'Color', [0.9, 0.9, 0.9]);

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10, 10, 4.3, 3.5]);
set(gcf, 'PaperSize', [4.3, 3.5]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [1, 0.8, 3, 2.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 6);
set(gca, 'XColor', [0.9, 0.9, 0.9]);
set(gca, 'YColor', [0.9, 0.9, 0.9]);


xlim ([-400, 400]);
set(gca, 'XTick', (-400:200:400));
set(gca, 'XMinorTick', 'off');
set(gca, 'XGrid', 'on');

ylim ([-0.5, 1]);
set (gca, 'YTick', (-0.5:0.25:1));
set (gca, 'YMinorTick', 'off');
set (gca, 'YGrid', 'on');
set(gca, 'GridLineStyle', '-.');

xl_hdl = xlabel ('$\tau$ / ns', 'Interpreter', 'LaTeX', 'FontSize', 8);
xl_pos = get (xl_hdl, 'Position');
set (xl_hdl, 'Position', xl_pos + [0, 0.1, 0]);

ylbl_hdl = ylabel ('$p_{0}(\tau) / \left\|p_{0}\right\|_{\infty}$', 'Units', 'centimeters', 'Interpreter', 'LaTeX', 'FontSize', 8);
ylbl_pos = get(ylbl_hdl, 'Position');
set (ylbl_hdl, 'Position', ylbl_pos + [0.55, 0, 0]);

print('-r600', '-dpdf', 'pressure_initial_gray.pdf');