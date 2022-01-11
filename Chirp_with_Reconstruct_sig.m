%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lizette Tovar, lizette.tovar-torres@uni-ulm.de %%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Example of chirp genereation example                                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 10.01.22 correction when calculating correlation of complex signal  %%%
%%% 15.12.21 Creating interferer function                               %%%
%%% 09.12.21 include interferer                                         %%%
%%% 08.12.21 starting version                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% include additive white gaussian noise
add_AWGN       = true;           % select: true or false 

% apply window to time data
include_window = true;           % select: true or false 

% chirps to evaluate
number_of_chips = 1;              % select: 1 or 128     

% signal type 
signal_type = 'complex';          % select: 'real' or 'complex'

% plot axis in selected scale
plot_axis = 'm_m/s';              % select: 'samples', 'Hz', 'm_m/s'
 
%% configuration parameters

cfg.c = physconst('LightSpeed');  % speed of light

% target parametes (in the example 2 targets are considered)
cfg.target.r   = [10 5  25];  % in m
cfg.target.v   = [10 20 30]; % in m/s
cfg.target.rcs = [10 5  2]; % in dbsm
cfg.target.phase = rand(1, length(cfg.target.r)) * 2 * pi;

% parametes of sensor 1 (victim)
cfg.radar.fc  = 79e9;  % in Hz
cfg.radar.fs  = 30e6;  % in Hz 
cfg.radar.B   = 1e9;   % in Hz
cfg.radar.Brx = 15e6;  % in Hz
cfg.radar.Tc  = 20e-6; % in s

cfg.radar.Tr  = 25e-6; % in s
cfg.radar.Tb  = 20e-3; % in s
cfg.radar.nr_chirp = number_of_chips;

% parametes of sensor 2 (interferer)
cfg.interf.fc = 79e9;   % in Hz
cfg.interf.B  = 1e9;    % in Hz
cfg.interf.Tc = 75e-6;  % in s
cfg.interf.Tr  = 80e-6; % in s
cfg.interf.Tb  = 32e-3; % in s
cfg.interf.nr_chirp = 1;
cfg.interf.r = [20 25];      % in m
cfg.interf.v = [10 15];      % in m/s
cfg.interf.shift = 0;2e-6;   % in s


%% select target to consider
i_tar = 1; 


%%  create time vector
Tc = ceil( cfg.radar.Tc * cfg.radar.fs );  % number of samples for each ramp
nr_chirp =  cfg.radar.nr_chirp;            % number of chirps
Tc_t = (0:Tc-1) / cfg.radar.fs;
time = repmat( ( Tc_t - Tc_t(round(Tc/2)) ), [nr_chirp,1]); %[nr_chirp X Tc]


%% create phase of signal due to target 

R_ramp =    cfg.target.r(i_tar) ...                               % is the range of the simulated target 
         + (0:nr_chirp-1)' * cfg.target.v(i_tar) * cfg.radar.Tr;  % is the distance the target travels 

% Note: remember thao is 2*R_ramp/c for this case

phi_targ = 2*pi*  (...
                    2 * cfg.radar.fc * R_ramp / cfg.c ...                                         % constant phase term 2*fc*R/C
                    + (2 * cfg.radar.fc * cfg.target.v(i_tar) / cfg.c ...                         % doppler term (time dependant) +(2*fc*v/C
                    + 2 * cfg.radar.B * R_ramp / (cfg.radar.Tc * cfg.c) ) ...                     % range term (time dependant)   + 2B*R/TC)
                    .* time + ...                                                                 % time dependency               *t
                    2 * cfg.radar.B * cfg.target.v(i_tar) / (cfg.radar.Tc * cfg.c) .* time.^2 ... % range doppler coupling        +2*B*v/TC *t^2
                    );

                
%% create amplitude of signal due to target      
cfg.radar.Pt =  10;           % in dB   Transmitted Power     
cfg.radar.Gt =  20;           % in dB   Gain of transmitting antenna
cfg.radar.Gr =  20;           % in dB   Gain of receiving antenna
cfg.radar.gain_receiver = 33; % in dB   Receivers gain due to link budget

% calulate Pr = Pt * Gt * (Gr*rec_gain) * lambda^2 / (4pi)^3 * rcs / R^4
rcs_scale = cfg.radar.Pt + cfg.radar.Gt + cfg.radar.Gr + cfg.target.rcs(i_tar) + 2*10*log10( cfg.c/cfg.radar.fc ) - 10*log10( (4*pi)^3 * cfg.target.r(i_tar).^4 ) + cfg.radar.gain_receiver; % [in dBW]
rcs_scale = 10.^(rcs_scale/20); % in V.  P ~ u^2 -> sqrt(P) ~ u

                
%% create time signal

if strcmp(signal_type, 'real')
    % -> real signal 
    m_time_targ = sqrt(2) * rcs_scale * cos( phi_targ + cfg.target.phase(i_tar));  
elseif strcmp(signal_type, 'complex')
    % -> complex signal
    m_time_targ = sqrt(2) * 1/2 * rcs_scale * exp( 1i * phi_targ ) * exp(1i * cfg.target.phase(i_tar)) ; 
else
    error('------>  wrong selection of "signal_type"')
end



%% add AWGN
if add_AWGN
    cfg.radar.NF = 15;              % in dB   Receives noise figure

    
    g_chain = cfg.radar.Gr + cfg.radar.gain_receiver;

    % calculate awgn noise power sigma^2 = K*Tout*Brx = K*(Te*G)*B= K*((Fchain-1)To*G)*B
    n_sigma = physconst('Boltzmann') * (10^(cfg.radar.NF/10)-1) * 290 * cfg.radar.Brx  * 10^(g_chain/10); % add awgn noise power in dBW according to Philipps radar
                                                                                                          % N = K*Tout*Brx = K*(Te*G)*B= K*((Fchain-1)To*G)*B
                                                                                                          % Fchain    -> receiver noise figure (including everything from antenna to adc) 
                                                                                                          % To=290K   -> room temperature
                                                                                                          % G         -> gain of the complete receiver chain
                                                                                                          % B         -> IF low pass filter cut off                                                                                       

    noise_AWGN = sqrt(n_sigma) * randn(size(m_time_targ));

    % add additive withe gaussian noise
     m_time_targ = m_time_targ + noise_AWGN;
end


%% apply window (optional)
if include_window
    m_time_targ = repmat( hann(size(m_time_targ,2))',[nr_chirp,1]) .* m_time_targ ; % apply window in sample   direction;
    m_time_targ = repmat( hann(size(m_time_targ,1)) ,[1,Tc])       .* m_time_targ ; % apply window in velocity direction;
end


%% plot signal due to target in time
num_fig = 10;
chirp_selec = ceil(cfg.radar.nr_chirp/2);    % if more than one chirp is configured and if the window is applied, then set chirp_selec ~= 1;


figure(num_fig);  clf(num_fig); hold on;
plot(time,real(m_time_targ(chirp_selec,:)))
xlabel('time')
ylabel('amplitude in V')
title(['time signal over chirp no. ' num2str(chirp_selec)])

% plot time matrix in 3D if multiple chirps are considered
if size(m_time_targ,1) > 1
    figure(num_fig+1);  clf(num_fig+1);
    surf(real(m_time_targ),'EdgeColor','none')
    axis tight; view([0,90])
    xlabel('samples')
    ylabel('chirps')
    title('time signal over all the chirps')
end


%% find target's range
num_samples = Tc;
fft_ran = fft(m_time_targ, num_samples, 2)/num_samples;       % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
fft_ran = 10*log10(abs( fft_ran ).^2);                        % in dB
        
f_r     = (0: num_samples -1) * cfg.radar.fs / (num_samples); % x axis in frequency
range   = (cfg.radar.Tc * cfg.c) / (2 * cfg.radar.B) * f_r;   % x axis in m
    
% select positive half of the spectrum because targest with negative range do not make sense  
fft_ran = fft_ran(chirp_selec, 1:round(num_samples/2));
f_r     = f_r(1:round(num_samples/2));
range   = range(1:round(num_samples/2));
fft_ran = fft_ran + 3; % add 3dB because just half of the spectrum is considered


figure(num_fig+2);  clf(num_fig+2)

if strcmp(plot_axis, 'samples')
    % -> in samples
    plot(fft_ran);        xlabel('samples');    ylabel('power in dB');  grid on; title(['fft signal over chirp no. ' num2str(chirp_selec)])
elseif strcmp(plot_axis, 'Hz')
    % -> in Hz
    plot(f_r,fft_ran);    xlabel('freq in Hz'); ylabel('power in dB');  grid on; title(['fft signal over chirp no. ' num2str(chirp_selec)])
elseif strcmp(plot_axis, 'm_m/s')
    % -> in meters
    plot(range,fft_ran);  xlabel('range in m'); ylabel('power in dB');  grid on; title(['fft signal over chirp no. ' num2str(chirp_selec)])
    hold on; plot(cfg.target.r(i_tar), 20*log10(rcs_scale),'xr') % expected target
else
    error('------>  wrong selection of "plot_axis"')
end


%% find the targets velocity if multiple chirps are evaluated
if size(m_time_targ,1) > 1
    
    % 1° FFT across samples
    fft2_ran = fft(m_time_targ, num_samples, 2)/num_samples;    % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
    fft2_ran = fft2_ran(:,(1:round(num_samples/2)));            % take the positive part of the spectrum 
    
    % 2° FFT across chirps
    num_ramps = cfg.radar.nr_chirp;
    fft2_ran_vel = fft(fft2_ran, num_ramps, 1) / num_ramps;       % fft_vel( , ,1) along the columns of time_data: fft across the chirps in the 1°, 2°,... bin or sample. Gives fv doppler freq
    fft2_ran_vel = fftshift(fft2_ran_vel,1);                      % spectrum centered at 0Hz
    
    
    fft2_ran_vel = 10*log10(abs( fft2_ran_vel ).^2);              % in dB
    fft2_ran_vel = fft2_ran_vel + 3;                              % add 3dB because just half of the range spectrum is considered
    
    
    f_v = (-num_ramps/2:  num_ramps/2-1) * (1/cfg.radar.Tr) / num_ramps; % y axis in frequency
    velocity =  cfg.c / ( 2 * cfg.radar.fc ) * f_v;                      % y axis in m/s
    
    
    figure(num_fig+3);  clf(num_fig+3);
    if strcmp(plot_axis, 'samples')
        % -> in chirps Vs samples
        surf(fft2_ran_vel,'EdgeColor','none');                   axis tight; view([0,90]);    xlabel('samples');    ylabel('chirps');           title('2fft signal over all the chirps')
    elseif strcmp(plot_axis, 'Hz')
        % -> in Hz
        surf(f_r, f_v, fft2_ran_vel,'EdgeColor','none');         axis tight; view([0,90]);    xlabel('fr Hz');      ylabel('fv Hz');            title('2fft signal over all the chirps') 
    elseif strcmp(plot_axis, 'm_m/s')
        % -> in m/s Vs m
        surf(range, velocity, fft2_ran_vel,'EdgeColor','none');  axis tight; view([0,90]);    xlabel('range in m'); ylabel('velocity in m/s');  title('2fft signal over all the chirps')  
        hold on; plot3(cfg.target.r(i_tar), cfg.target.v(i_tar), 20*log10(rcs_scale),'xr') % expected target
    else
        error('------>  wrong selection of "plot_axis"')
    end
   
end

%% Interferer introduction

interferer_1= createInterferer(cfg,time, signal_type, 1);
interferer_2= createInterferer(cfg,time, signal_type, 2);

rx_sig_1 = m_time_targ + interferer_1;
rx_sig_2 = interferer_2;

step=0.01;
max_realpeak = 0;
max_imagpeak = 0; 

for n = (1:step:360)*pi/180
   processed_real_sig = interferer_2 .* exp( 1i * n );
   realcomp_data_peak = max(xcorr(real(rx_sig_1), real(processed_real_sig)));
   
   %save phase 'n' for which we get maximum correlation of the real part
   if realcomp_data_peak > max_realpeak 
        max_realpeak = abs(realcomp_data_peak);
        max_realpeak_phase = n;
   end
end

for m = (1:step:360)*pi/180
   processed_imag_sig = interferer_2 .* exp( 1i * m );
   imagcomp_data_peak = max(xcorr(imag(rx_sig_1), imag(processed_imag_sig)));

    %save phase 'm' for which we get maximum correlation of the imag part
   if imagcomp_data_peak > max_imagpeak 
        max_imagpeak  = angle(imagcomp_data_peak);
        max_imagpeak_phase = m;
   end
end

final_processed_realsig = interferer_2 .* exp( 1i * max_realpeak_phase );
final_processed_imagsig = interferer_2 .* exp( 1i * max_imagpeak_phase );
reconstruct_sig = rx_sig_1 - final_processed_realsig - final_processed_imagsig;

figure; 
subplot(2,1,1); hold on; 
plot(real(m_time_targ))                                          % clean signal
plot(real(rx_sig_1))                                             % signal with interference
plot(real(reconstruct_sig))                                      % signal after mitigation
title('real part'); legend('T', 'T+I', 'mitigated')
subplot(2,1,2); hold on;
plot(imag(m_time_targ))                                          % clean signal
plot(imag(rx_sig_1))                                             % signal with interference
plot(imag(reconstruct_sig))                                      % signal after mitigation
title('imaginary part'); legend('T', 'T+I', 'mitigated')

figure; 
subplot(2,1,1); hold on;
plot(abs(m_time_targ))                                          % clean signal
plot(abs(rx_sig_1))                                             % signal with interference
plot(abs(reconstruct_sig))                                      % signal after mitigation
title('magnitude'); legend('T', 'T+I', 'mitigated')
subplot(2,1,2); hold on;                                            
plot(unwrap(angle(m_time_targ)))                                          % clean signal
plot(unwrap(angle(rx_sig_1)))                                             % signal with interference
plot(unwrap(angle(reconstruct_sig)))                                      % signal after mitigation
title('phase'); legend('T', 'T+I', 'mitigated')


%we have our phase 'max_peak_phase' which when multiplied with 'rx_sig_2' nullifies interference
% final_processed_sig = interferer_2 .* exp( 1i * max_peak_phase_2 );
% comp_data_peak = abs(xcorr(rx_sig_1, final_processed_sig));
% figure("Name", "Maximum Correlation with new found phase")
% plot(comp_data_peak)

% final_rx_sig(final_rx_sig==0)=1;

% reconstruct_sig = rx_sig_1 .* conj(final_rx_sig);
figure("Name", "Reconstructed Signal")
plot(real(reconstruct_sig))
% figure("Name", "RX Signal 1")
% plot(real(rx_sig_1))
% figure("Name", "Final RX Signal")
% plot(real(final_rx_sig))



% %%-> plot interferer in time
% figure(num_fig);  hold on;
% plot(time,real(interferer_1),'r')       % interferer
% xline(e,':r')                         % interferer center
% legend('target','interferer','interf center')
% 
% 
% % -> fft of interferer
% fft_ran_i = fft(interferer_1, num_samples, 2)/num_samples;       % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
% fft_ran_i = 10*log10(abs( fft_ran_i ).^2);                     % in dB
% 
% % -> fft of interferer + target
% fft_ran_i_t = fft(rx_sig_1, num_samples, 2)/num_samples;  % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
% fft_ran_i_t = 10*log10(abs( fft_ran_i_t ).^2);                            % in dB
%     
% % select positive half of the spectrum because targest with negative range do not make sense  
% fft_ran_i = fft_ran_i(chirp_selec, 1:round(num_samples/2));
% fft_ran_i = fft_ran_i + 3; % add 3dB because just half of the spectrum is considered
% fft_ran_i_t = fft_ran_i_t(chirp_selec, 1:round(num_samples/2));
% fft_ran_i_t = fft_ran_i_t + 3; % add 3dB because just half of the spectrum is considered
% 
% figure(num_fig+2); hold on
% 
% if strcmp(plot_axis, 'samples')
%     % -> in samples
%     plot(fft_ran_i);     
%     plot(fft_ran_i_t);  
% elseif strcmp(plot_axis, 'Hz')
%     % -> in Hz
%     plot(f_r,fft_ran_i);
%     plot(f_r,fft_ran_i_t); 
% elseif strcmp(plot_axis, 'm_m/s')
%     % -> in meters
%     plot(range,fft_ran_i);
%     plot(range,fft_ran_i_t);
% else
%     error('------>  wrong selection of "plot_axis"')
% end
% legend('target','expected peak', 'interferer','tar+interf')


%% find target's range with reconstruct signal
num_samples = Tc;
fft_ran = fft(reconstruct_sig, num_samples, 2)/num_samples;       % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
fft_ran = 10*log10(abs( fft_ran ).^2);                        % in dB
        
f_r     = (0: num_samples -1) * cfg.radar.fs / (num_samples); % x axis in frequency
range   = (cfg.radar.Tc * cfg.c) / (2 * cfg.radar.B) * f_r;   % x axis in m
    
% select positive half of the spectrum because targest with negative range do not make sense  
fft_ran = fft_ran(chirp_selec, 1:round(num_samples/2));
f_r     = f_r(1:round(num_samples/2));
range   = range(1:round(num_samples/2));
fft_ran = fft_ran + 3; % add 3dB because just half of the spectrum is considered


figure(num_fig+7);  clf(num_fig+7)

if strcmp(plot_axis, 'samples')
    % -> in samples
    plot(fft_ran);        xlabel('samples');    ylabel('power in dB');  grid on; title(['fft reconstruct signal over chirp no. ' num2str(chirp_selec)])
elseif strcmp(plot_axis, 'Hz')
    % -> in Hz
    plot(f_r,fft_ran);    xlabel('freq in Hz'); ylabel('power in dB');  grid on; title(['fft reconstruct signal over chirp no. ' num2str(chirp_selec)])
elseif strcmp(plot_axis, 'm_m/s')
    % -> in meters
    plot(range,fft_ran);  xlabel('range in m'); ylabel('power in dB');  grid on; title(['fft reconstruct signal over chirp no. ' num2str(chirp_selec)])
    hold on; plot(cfg.target.r(i_tar), 20*log10(rcs_scale),'xr') % expected target
else
    error('------>  wrong selection of "plot_axis"')
end





