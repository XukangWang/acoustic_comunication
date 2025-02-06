%% BFSK Communication System Simulation 
% Clear workspace and command window
clear; clc; close all;

%% Parameter Settings
fs = 192e3;         % Sampling rate 192 kHz
T_sym = 1e-3;       % Symbol duration 1 ms
t_sym = 0:1/fs:T_sym-1/fs;  % Time vector for a single symbol

f0 = 20e3;          % Carrier frequency for bit 0: 20 kHz
f1 = 40e3;          % Carrier frequency for bit 1: 40 kHz

numSymbols = 50;   % Number of symbols

%% Generate Random Bit Sequence
data = randi([0, 1], 1, numSymbols);

%% BFSK Modulation
% Generate a sine wave for each bit corresponding to its frequency and concatenate all symbols
tx_signal = [];
for k = 1:numSymbols
    if data(k) == 0
        sig = cos(2*pi*f0*t_sym);
    else
        sig = cos(2*pi*f1*t_sym);
    end
    tx_signal = [tx_signal, sig];
end

%% Underwater Acoustic Channel Simulation
% The underwater acoustic channel typically exhibits multipath effects. Here, we use a simple multipath model:
%  - First direct path, gain = 1, no delay
%  - Second path, delay of 2~3 samples, gain ≈ 0.6
%  - Third path, delay of 7~8 samples, gain ≈ 0.3
% Construct a discrete impulse response (length can be adjusted as needed)
h = [1, zeros(1,2), 0.6, zeros(1,4), 0.3];

% Simulate multipath effects using a filter
channel_signal = filter(h, 1, tx_signal);

%% Additive Noise (AWGN)
% Set signal-to-noise ratio (SNR) to 15 dB (modifiable as needed)
SNR_dB = 15;
rx_signal = awgn(channel_signal, SNR_dB, 'measured');

%% BFSK Demodulation (FFT Spectrum Analysis)
% Segment the received signal based on the number of samples per symbol
N_sym = length(t_sym);   % Number of samples per symbol
rx_matrix = reshape(rx_signal, N_sym, []);

% FFT frequency resolution is fs/N_sym = 192e3/192 = 1 kHz,
% Therefore:
%   - FFT bin index for 20 kHz: index_f0 = 20 + 1 = 21
%   - FFT bin index for 40 kHz: index_f1 = 40 + 1 = 41
index_f0 = round(f0 * N_sym / fs) + 1;  
index_f1 = round(f1 * N_sym / fs) + 1;  

detected = zeros(1, numSymbols);  % Initialize detection results

for k = 1:numSymbols
    symbol_signal = rx_matrix(:, k);   % Extract the k-th symbol signal
    
    % Compute FFT
    X = fft(symbol_signal);
    magX = abs(X);   % Compute magnitude spectrum
    
    % Decision based on single bin energy
    amp0 = magX(index_f0);
    amp1 = magX(index_f1);
    
    % If 20 kHz component is stronger, detect as bit 0; otherwise, detect as bit 1
    if amp0 > amp1
        detected(k) = 0;
    else
        detected(k) = 1;
    end
end

%% Compute Bit Error Rate (BER)
numErrors = sum(data ~= detected);
BER = numErrors / numSymbols;
fprintf('Bit Error Rate (BER) = %f\n', BER);

%% Plot Results
figure;
subplot(4,1,1);
plot(tx_signal);
title('Transmitted BFSK Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(4,1,2);
plot(channel_signal);
title('Signal After Underwater Acoustic Channel');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(4,1,3);
plot(rx_signal);
title('Received Signal with Added Noise');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(4,1,4);
stem(data, 'filled'); hold on;
stem(detected, 'r', 'filled');
title('Original Bits (Blue) vs Detected Bits (Red)');
xlabel('Symbol Index');
ylabel('Bit Value');
legend('Original Bits','Detected Bits');
