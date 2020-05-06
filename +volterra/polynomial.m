function pressure_output = polynomial( pressure_input, f_s, delta_z, a, b, order_N )
%
% author: Martin Schiffner
% date: 2009-04-14
% modified: 2020-05-06

% internal parameters
N_input = numel(pressure_input);
atten_td = 150;       % attenuation in time-domain at filter-length (dB)

%create the Gauss filter's transfer function
%choose length to be an odd number, linear phase
M_filter_gauss = ceil(f_s * sqrt(atten_td * abs(a) * delta_z * log(10) / 5));
N_filter_gauss = 2 * M_filter_gauss + 1;
N_fft = N_filter_gauss + N_input - 1;     %take length of convolution into account

nyquist_index = ceil(N_fft / 2);
Omega_filter = 2 * pi * ((nyquist_index - N_fft):(nyquist_index - 1)) / N_fft;
filter_gauss_fft = ifftshift(exp(-abs(a) * (Omega_filter * f_s).^2 * delta_z)) .* exp(-1j * 2 * pi * (0:(N_fft - 1)) * M_filter_gauss / N_fft);

%sample Gauss filter more densly for filtering the periodically extended signal
%same filter length, longer input signal, linear phase
N_fft_2 = N_filter_gauss + 2 * N_input - 1;
nyquist_index_2 = ceil(N_fft_2 / 2);
Omega_filter_2 = 2 * pi * ((nyquist_index_2 - N_fft_2):(nyquist_index_2 - 1)) / N_fft_2;
filter_gauss_fft_2 = ifftshift(exp(-abs(a) * (Omega_filter_2 * f_s).^2 * delta_z)) .* exp(-1j * 2 * pi * (0:(N_fft_2 - 1)) * M_filter_gauss / N_fft_2);

%allocate memory for output signals
signals_output = cell(1, order_N);

%compute output of linear system
signals_output{1,1} = ifft(filter_gauss_fft .* fft(pressure_input, N_fft), 'symmetric');
pressure_output = signals_output{1,1};
    
%compute outputs of homogeneous systems of order >= 2
for l = 2:order_N
        
    signals_output{1,l} = volterra.Hn_structure_IIR_6_simplified( pressure_input, filter_gauss_fft, filter_gauss_fft_2, f_s, a, b, signals_output(1, 1:(l - 1)) );
    pressure_output = pressure_output + signals_output{1,l};
end

%cut-off output signal
pressure_output = pressure_output((M_filter_gauss + 1):(end - M_filter_gauss));

%apply window function
pressure_output = tukeywin(N_input, 0.1)' .* pressure_output;

end % function pressure_output = polynomial( pressure_input, f_s, delta_z, a, b, order_N )
