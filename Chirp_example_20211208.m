%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lizette Tovar, lizette.tovar-torres@uni-ulm.de %%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Example of chirp genereation example                                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 15.12.21 Creating interferer function & modifying 'g'               %%%
%%% 09.12.21 include interferer                                         %%%
%%% 08.12.21 starting version                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% include additive white gaussian noise
add_AWGN       = true;           % select: true or false 

% apply window to time data
include_window = false;           % select: true or false 

% chirps to evaluate
number_of_chips = 1;              % select: 1 or 128     

% signal type 
signal_type = 'real';          % select: 'real' or 'complex'

% plot axis in selected scale
plot_axis = 'm_m/s';              % select: 'samples', 'Hz', 'm_m/s'
 
%% configuration parameters

cfg.c = physconst('LightSpeed');  % speed of light

% target parametes (in the example 2 targets are considered)
cfg.target.r   = [10  5];  % in m
cfg.target.v   = [10 10]; % in m/s
cfg.target.rcs = [10 5]; % in dbsm
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
cfg.interf.r = 10;      % in m
cfg.interf.v = 10;      % in m/s
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

                
%% create amplitud of signal due to target      
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


%% apply findow (optional)
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



%% create interferer

int_chirp_idx =  1;   % number of affecting interfering chirp

int_slp = cfg.interf.B/cfg.interf.Tc; % interferer slope
radar_slp = cfg.radar.B/cfg.radar.Tc; % victim slope

%-> samples of interference duration
int_duration = floor(2*cfg.radar.Brx / abs(int_slp-radar_slp) *cfg.radar.fs);


% -> find terms related to the phase
R_int= cfg.interf.r + ...                                % is the initial distance within interf and victim
       (int_chirp_idx-1) * cfg.interf.v * cfg.interf.Tr; % is the distance the interf travels each interferer ramp 

thao_int = cfg.interf.shift + R_int/cfg.c;        

b = cfg.interf.fc -  cfg.radar.fc - cfg.interf.fc*cfg.interf.v/cfg.c - int_slp*thao_int + int_slp*thao_int*cfg.interf.v/cfg.c; % b = fci -fc -fci*vi/c -ui*ri/c +ui*ri/c*vi/c
a = int_slp/2 - radar_slp/2 - int_slp*cfg.interf.v/cfg.c;                                                                      % a = ui/2 -u/2 -ui*vi/c +0
c = -cfg.interf.fc*thao_int + int_slp*(thao_int).^2/2;                                                                         % c = -fci*ri/c +ui/2*(ri/c)^2 
e = -b/(2*a);



% -> get frequency ramps of victim and interferer  ft= fc + B/Tc*t
time = time(1,:);  % consider just one chirp

time_int = circshift(time,round(thao_int * cfg.radar.fs)-1);
time_int(1:round(thao_int * cfg.radar.fs)) = NaN;
time_int(end + round(thao_int * cfg.radar.fs):end) = NaN;


f_t   = cfg.radar.fc + cfg.radar.B/cfg.radar.Tc * time;        
f_int = cfg.interf.fc + cfg.interf.B/cfg.interf.Tc * time_int;
figure(num_fig+4);  clf(num_fig+4); hold on;           
plot(time,f_t,'k')                    % victim ramp
plot(time,f_t + cfg.radar.Brx,'--k')  % upper receiver bandwidth 
plot(time,f_t - cfg.radar.Brx,'--k')  % lower receiver bandwidth 
plot(time,f_int, 'r')               % interferer ramp
xline(e,':r')                         % interferer center
xlabel('time'); ylabel('freq in Hz')
legend('victim','up Brx','low Brx','interferer','interf center')
grid on


% -> get sample where the interferer center interferer center is placed
tmp1 = f_t - f_int;
tmp2 = tmp1 > 0;
int_center = find((diff(tmp2)==1 | diff(tmp2)==-1)); %index of interferer intersections
int_center(isnan(tmp1(int_center+1)))=[];            %remove the indexes related to NaNs (+ NaN)
int_center(isnan(tmp1(int_center)))=[];              %remove the indexes related to NaNs (NaN +)
% int_center = int_center;

if isempty(int_center)
    error('---> set shorter cfg.interf.shift')
end

% -> time vector for interferer
int_samples = (-floor(int_duration/2):1:floor(int_duration/2)) + int_center; 
time_int = NaN(size(time));
time_int(int_samples) = time(int_samples);



% -> calculate phase of interferer phase_int = 2pi*( b*t + a*t^2 + c )
% phase_int = 2*pi* (b.*time_int + a*time_int.^2 + c);
phase_int = 2*pi* (b.*time_int + a*time_int.^2 + c);


%-> create amplitud of signal due to interferer      
cfg.interf.Gt =  20;           % in dB   Gain of transmitting antenna

% calulate Pr = Pt * Gt * (Gr*rec_gain) * lambda^2 / (4pi * R)^2
interf_lev = cfg.radar.Pt +  cfg.interf.Gt + cfg.radar.Gr + 10*log10(( cfg.c/cfg.radar.fc./(4*pi*cfg.interf.r) ).^2) + cfg.radar.gain_receiver; % [in dBW];
interf_lev = 10.^(interf_lev/20); % in V.  P ~ u^2 -> sqrt(P) ~ u



% -> create time signal
if strcmp(signal_type, 'real')
    % -> real signal 
    m_time_int = sqrt(2) * interf_lev * cos( phase_int);
elseif strcmp(signal_type, 'complex')
    % -> complex signal
   m_time_int = sqrt(2) * 1/2 * interf_lev * exp( 1i * phase_int ) ; 
else
    error('------>  wrong selection of "signal_type"')
end
m_time_int(isnan(m_time_int)) = 0;



%-> plot interferer in time
figure(num_fig);  hold on;
plot(time,real(m_time_int),'r')       % interferer
xline(e,':r')                         % interferer center
legend('target','interferer','interf center')



% -> fft of interferer
fft_ran_i = fft(m_time_int, num_samples, 2)/num_samples;       % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
fft_ran_i = 10*log10(abs( fft_ran_i ).^2);                     % in dB

% -> fft of interferer + target
fft_ran_i_t = fft(m_time_int + m_time_targ, num_samples, 2)/num_samples;  % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
fft_ran_i_t = 10*log10(abs( fft_ran_i_t ).^2);                            % in dB
    
% select positive half of the spectrum because targest with negative range do not make sense  
fft_ran_i = fft_ran_i(chirp_selec, 1:round(num_samples/2));
fft_ran_i = fft_ran_i + 3; % add 3dB because just half of the spectrum is considered
fft_ran_i_t = fft_ran_i_t(chirp_selec, 1:round(num_samples/2));
fft_ran_i_t = fft_ran_i_t + 3; % add 3dB because just half of the spectrum is considered

figure(num_fig+2); hold on

if strcmp(plot_axis, 'samples')
    % -> in samples
    plot(fft_ran_i);     
    plot(fft_ran_i_t);  
elseif strcmp(plot_axis, 'Hz')
    % -> in Hz
    plot(f_r,fft_ran_i);
    plot(f_r,fft_ran_i_t); 
elseif strcmp(plot_axis, 'm_m/s')
    % -> in meters
    plot(range,fft_ran_i);
    plot(range,fft_ran_i_t);
else
    error('------>  wrong selection of "plot_axis"')
end
legend('target','expected peak', 'interferer','tar+interf')








