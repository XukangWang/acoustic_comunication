%% BFSK 通信系统仿真
% 清空工作区和命令窗口
clear; clc; close all;

%% 参数设置
fs = 192e3;         % 采样率 192 kHz
T_sym = 1e-3;       % 每个符号时长 1 ms
t_sym = 0:1/fs:T_sym-1/fs;  % 单个符号的时间向量

f0 = 20e3;          % 比特 0 对应的载波频率 20 kHz
f1 = 40e3;          % 比特 1 对应的载波频率 40 kHz

numSymbols = 50;   % 符号数

%% 生成随机比特序列
data = randi([0, 1], 1, numSymbols);

%% BFSK 调制
% 对每个比特，生成对应频率的正弦波，并将所有符号级联
tx_signal = [];
for k = 1:numSymbols
    if data(k) == 0
        sig = cos(2*pi*f0*t_sym);
    else
        sig = cos(2*pi*f1*t_sym);
    end
    tx_signal = [tx_signal, sig];
end

%% 水声信道仿真
% 水声信道一般具有多径效应，这里采用一个简单的多径模型
% 例如，假设信道的脉冲响应包含三个路径：
%  - 第一条直达路径，增益 = 1，无延时
%  - 第二条路径，延时 2~3 个采样点，增益约为 0.6
%  - 第三条路径，延时 7~8 个采样点，增益约为 0.3
% 构造一个离散脉冲响应（单位冲激响应长度可以根据实际情况调整）
h = [1, zeros(1,2), 0.6, zeros(1,4), 0.3];

% 通过滤波器模拟多径效应
channel_signal = filter(h, 1, tx_signal);

%% 加性噪声（AWGN）
% 设定信噪比（SNR）为 20 dB（可根据需要修改）
SNR_dB = 20;
rx_signal = awgn(channel_signal, SNR_dB, 'measured');

%% BFSK 解调（FFT频谱分析）
% 将接收信号按照每个符号的采样点数分段
N_sym = length(t_sym);   % 每个符号的采样点数
rx_matrix = reshape(rx_signal, N_sym, []);

% FFT 频率分辨率为 fs/N_sym = 192e3/192 = 1 kHz，
% 因此：
%   - 20 kHz 对应的 FFT bin 索引：index_f0 = 20 + 1 = 21
%   - 40 kHz 对应的 FFT bin 索引：index_f1 = 40 + 1 = 41
index_f0 = round(f0 * N_sym / fs) + 1;  
index_f1 = round(f1 * N_sym / fs) + 1;  

detected = zeros(1, numSymbols);  % 初始化检测结果

for k = 1:numSymbols
    symbol_signal = rx_matrix(:, k);   % 取出第 k 个符号信号
    
    % 计算 FFT
    X = fft(symbol_signal);
    magX = abs(X);   % 求幅值谱
    
    % 采用单个 bin 的能量进行判决
    amp0 = magX(index_f0);
    amp1 = magX(index_f1);
    
    % 如果 20 kHz 分量能量更强，判决为比特 0，否则为比特 1
    if amp0 > amp1
        detected(k) = 0;
    else
        detected(k) = 1;
    end
end

%% 计算误码率（BER）
numErrors = sum(data ~= detected);
BER = numErrors / numSymbols;
fprintf('误码率 = %f\n', BER);

%% 绘制结果
figure;
subplot(4,1,1);
plot(tx_signal);
title('发送的 BFSK 信号');
xlabel('采样点');
ylabel('幅值');

subplot(4,1,2);
plot(channel_signal);
title('经过水声信道后的信号');
xlabel('采样点');
ylabel('幅值');

subplot(4,1,3);
plot(rx_signal);
title('加噪声后的接收信号');
xlabel('采样点');
ylabel('幅值');

subplot(4,1,4);
stem(data, 'filled'); hold on;
stem(detected, 'r', 'filled');
title('原始比特（蓝色）与检测比特（红色）');
xlabel('符号序号');
ylabel('比特值');
legend('原始比特','检测比特');
