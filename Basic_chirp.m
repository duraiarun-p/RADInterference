%% Initialization
clc; clear;

%% Chirp sequence parameters
chirp_B     = 1e9;                   % Sweep bandwidth
chirp_fc    = 77.5e9;                % Sweep starting frequency
chirp_T     = 20e-6;                 % Time period of one chirp or sweep time
chirp_L     = 1;                     % Chirp cycles
chirp_slope = chirp_B/chirp_T;       % f=at; Slope of chirp
c           = 3e8;                   % Speed of light
R           = 10;                    % Target distance
v           = 0;                     % Velocity of the target

%% ADC parameters
sf = chirp_B;                       % sampling frequency
bin_size = sf/10e6;                 % Bin size of each sample; Sample intervals

%% Time and frequency axis created for plots
t_axis = (0: 1/sf :chirp_L * chirp_T - 1/sf);
% t        = 0 : 1/sf : chirp_T;    % Time vector

%% Chirp generation
chirp_signal = exp(1i * pi * chirp_slope * t_axis.^2);
% freq_T = (chirp_fc + (chirp_slope*t));
% chirp_signal = cos(2*pi * freq_T.*t);

chirp_signal_buckets = reshape(chirp_signal(1:floor(length(chirp_signal)/bin_size)*bin_size), [bin_size, floor(length(chirp_signal)/bin_size)]);

f_axis = linspace(chirp_fc - chirp_B/2, chirp_fc + chirp_B/2, bin_size);
t_axis = linspace(0, chirp_T * chirp_L, floor(length(chirp_signal)/bin_size));

figure(1)
imagesc(t_axis, f_axis, 20*log10(abs(fft(chirp_signal_buckets))));
set(gca,'YDir','normal')
colorbar
axis tight

%% Time Delay for received signal
Td = (2/c)*(R+v*t_axis);

%% Received Chirp Signal

% r_chirp_signal = exp(1i * pi * chirp_slope * t_axis.^2);
% r_chirp_T_buckets = chirp_T - Td;
% 
% r_chirp_signal_buckets = reshape(r_chirp_signal(1:floor(length(r_chirp_signal)/bin_size)*bin_size), [bin_size, floor(length(r_chirp_signal)/bin_size)]);
% 
% f_axis = linspace(chirp_fc - chirp_B/2, chirp_fc + chirp_B/2, bin_size);
% t_axis = linspace(0, chirp_T * chirp_L, floor(length(r_chirp_signal)/bin_size));
% 
% figure(2)
% imagesc(t_axis, f_axis, 20*log10(abs(fft(r_chirp_signal_buckets))));
% set(gca,'YDir','normal')
% colorbar
% axis tight
