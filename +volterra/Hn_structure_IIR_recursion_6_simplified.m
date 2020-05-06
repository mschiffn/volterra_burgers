function y_out = Hn_structure_IIR_recursion_6_simplified(input_signals, filter_fft, filter_fft_2, f_s, a, b, y_out_prec)
%compute output signal for homogeneous Volterra-System of order
%order_N
%input_signals: cell array with 2 input signals: the first input signal
%               represents the first order_N - 1 input signals to a multilinear system
%               the second input signal is the order_N-th input signal to a
%               multilinear system of order order_N

%author: Martin Schiffner
%date: 2008-09-28

% general internal parameters
order_N = numel(y_out_prec) + 1;

%calculate intermediate signals
N_input = numel(input_signals{1,1});
N_fft = numel(filter_fft);
N_fft_2 = numel(filter_fft_2);

temp1 = filter(1 / f_s, [1, -1], input_signals{1,2});   %discrete integration using IIR filter (shift by + 1/2 sample!)

%interpolate both input signals by factor 2 before multiplication to avoid
%aliasing effects. Also compensate phase shift by + 1/2 sample!

%interpolate input1 manually
input1_fft = fft(input_signals{1,1});
nyquist_index = ceil(N_input / 2);
input1_interpolated = ifft(2 * [input1_fft(1:nyquist_index), zeros(1, N_input), input1_fft((nyquist_index + 1):end)], 'symmetric');

%interpolate temp1 manually and compensate for phase shift by +1/2 sample
%extend temp1 symmetrically to avoid edges in time-domain
temp1_padded_extended = [0, temp1, temp1(end:(-1):1), 0];   
temp1_padded_extended_fft = fft(temp1_padded_extended);
nyquist_index = N_input + 1;
temp1_padded_extended_interpolated_fft = 2 * [temp1_padded_extended_fft(1:nyquist_index), zeros(1, 2 * (N_input + 1)), temp1_padded_extended_fft((nyquist_index + 1):end)];
temp1_padded_extended_interpolated = ifft(temp1_padded_extended_interpolated_fft, 'symmetric');
temp1_padded_interpolated = temp1_padded_extended_interpolated(1:(2 * (N_input + 1)));
temp1_interpolated = temp1_padded_interpolated(2:(end - 1));

%multiplication @ 2*f_s to avoid aliasing
temp2_interpolated = temp1_interpolated .* input1_interpolated;

%use decimation to reduce sampling rate and to avoid aliasing
temp2_interpolated_fft = fft(temp2_interpolated);

%cut-off dft in frequency domain to filter out aliasing effects
nyquist_index = ceil(N_input / 2);
temp2_fft = [temp2_interpolated_fft(1:nyquist_index),temp2_interpolated_fft((N_input + nyquist_index + 1): end)] / 2;
temp2 = ifft(temp2_fft, 'symmetric');

if order_N == 2
    
    %calculate output of homogeneous first order system directly
    %use zero padding to avoid time-domain aliasing (N_fft is computed automatically)
    out1 = ifft(filter_fft .* fft(temp2, N_fft), 'symmetric');
        
elseif order_N > 2
    
    %apply window function to avoid edges
    temp2 = tukeywin(N_input, 0.1)' .* temp2;
    
    %use recurrence relation to calculate output of inhomogeneous system
    out1 = volterra.Hn_structure_IIR_recursion_6_simplified({input_signals{1,1}, temp2}, filter_fft, filter_fft_2, f_s, a, b, y_out_prec(1,1:(end - 1))); 
end

%use existing signal as output of homogeneous system
out2 = y_out_prec{1,end};
    
%apply Gauss filter to temp1_interpolated
%periodic extension not necessarily due to lowpass function of filter
out3_fft = filter_fft_2 .* fft([temp1_interpolated(1:2:end), temp1_interpolated(end:-2:1)], N_fft_2);

%interpolate output signals by factor 2 before multiplication to avoid aliasing
out2_interpolated = interpft(out2, 2 * numel(out2));

nyquist_index = ceil(N_fft_2 / 2);
out3_interpolated = ifft(2 * [out3_fft(1:nyquist_index), zeros(1, N_fft_2), out3_fft((nyquist_index + 1):end)], 'symmetric');

out23_interpolated = out2_interpolated .* out3_interpolated(1:(2 * N_fft));

%use decimation to reduce sampling rate and avoid aliasing
out23_interpolated_fft = fft(out23_interpolated);
out23_interpolated_fft = [out23_interpolated_fft(1:nyquist_index),out23_interpolated_fft((N_fft + nyquist_index + 1):end)] / 2;
out23 = ifft(out23_interpolated_fft, 'symmetric');

%compute output signal
y_out = b * (out1 - out23) / (2 * a);

%debugging:
%temp1_fft = fft(temp1, 2*numel(temp1));
%temp2_fft = fft(temp2, 2*numel(temp2));
%out1_fft = fft(out1, 2*numel(out1));
%out23_fft = fft(out23, 2*numel(out23));
%y_out_fft = fft(y_out, 2*numel(y_out));

%figure(1);
%plot(f_s * (0:(numel(temp1_fft) - 1)) / numel(temp1_fft), 20*log10(abs(temp1_fft) / max(abs(temp1_fft))));
%xlim([0 f_s/2]);
%figure(2);
%plot(2*f_s * (0:(numel(temp2_interpolated_fft) - 1)) / numel(temp2_interpolated_fft), 20*log10(abs(ifftshift(temp2_interpolated_fft)) / max(abs(temp2_interpolated_fft))));
%xlim([0 f_s]);
%figure(3);
%plot(f_s * (0:(numel(temp2_fft) - 1)) / numel(temp2_fft), 20*log10(abs(temp2_fft) / max(abs(temp2_fft))));
%xlim([0 f_s/2]);
%figure(4);
%plot(f_s * (0:(numel(out1_fft) - 1)) / numel(out1_fft), 20*log10(abs(out1_fft) / max(abs(out1_fft))));
%xlim([0 f_s / 2]);
%figure(5);
%plot(f_s * (0:(numel(out23_fft) - 1)) / numel(out23_fft), 20*log10(abs(out23_fft) / max(abs(out23_fft))));
%xlim([0 f_s/2]);
%figure(6);
%plot(f_s * (0:(numel(y_out_fft) - 1)) / numel(y_out_fft), 20*log10(abs(y_out_fft) / max(abs(y_out_fft))));
%xlim([0 f_s/2]);