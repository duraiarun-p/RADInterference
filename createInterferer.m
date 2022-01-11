function [interferer] = createInterferer(cfg, time, signal_type, select_int)

   int_chirp_idx =  1;   % number of affecting interfering chirp

int_slp = cfg.interf.B/cfg.interf.Tc; % interferer slope
radar_slp = cfg.radar.B/cfg.radar.Tc; % victim slope

%-> samples of interference duration
int_duration = floor(2*cfg.radar.Brx / abs(int_slp-radar_slp) *cfg.radar.fs);


% -> find terms related to the phase
R_int= cfg.interf.r(select_int) + ...                                % is the initial distance within interf and victim
       (int_chirp_idx-1) * cfg.interf.v(select_int) * cfg.interf.Tr; % is the distance the interf travels each interferer ramp 

thao_int = cfg.interf.shift + R_int/cfg.c;        

b = cfg.interf.fc -  cfg.radar.fc - cfg.interf.fc*cfg.interf.v(select_int)/cfg.c - int_slp*thao_int + int_slp*thao_int*cfg.interf.v(select_int)/cfg.c; % b = fci -fc -fci*vi/c -ui*ri/c +ui*ri/c*vi/c
a = int_slp/2 - radar_slp/2 - int_slp*cfg.interf.v(select_int)/cfg.c;                                                                      % a = ui/2 -u/2 -ui*vi/c +0
c = -cfg.interf.fc*thao_int + int_slp*(thao_int).^2/2;                                                                         % c = -fci*ri/c +ui/2*(ri/c)^2 
e = -b/(2*a);


% -> get frequency ramps of victim and interferer  ft= fc + B/Tc*t
time = time(1,:);  % consider just one chirp

time_int = circshift(time,round(thao_int * cfg.radar.fs)-1);
time_int(1:round(thao_int * cfg.radar.fs)) = NaN;
time_int(end + round(thao_int * cfg.radar.fs):end) = NaN;

%% Plot Interference signal

% num_fig = 20;
f_t   = cfg.radar.fc + cfg.radar.B/cfg.radar.Tc * time;        
f_int = cfg.interf.fc + cfg.interf.B/cfg.interf.Tc * time_int;
% figure(num_fig+4);  clf(num_fig+4); hold on;           
% plot(time,f_t,'k')                    % victim ramp
% plot(time,f_t + cfg.radar.Brx,'--k')  % upper receiver bandwidth 
% plot(time,f_t - cfg.radar.Brx,'--k')  % lower receiver bandwidth 
% plot(time,f_int, 'r')                 % interferer ramp
% xline(e,':r')                         % interferer center
% xlabel('time'); ylabel('freq in Hz')
% legend('victim','up Brx','low Brx','interferer','interf center')
% grid on


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
phase_int = 2*pi* (b.*time_int + a*time_int.^2 + c);



%-> create amplitud of signal due to interferer      
cfg.interf.Gt =  20;           % in dB   Gain of transmitting antenna

% calulate Pr = Pt * Gt * (Gr*rec_gain) * lambda^2 / (4pi * R)^2
interf_lev = cfg.radar.Pt +  cfg.interf.Gt + cfg.radar.Gr + 10*log10(( cfg.c/cfg.radar.fc./(4*pi*cfg.interf.r(select_int)) ).^2) + cfg.radar.gain_receiver; % [in dBW];
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

%% Generate interferer
interferer = m_time_int;

% 
% %%-> plot interferer in time
% figure(num_fig);  hold on;
% plot(time,real(m_time_int),'r')       % interferer
% xline(e,':r')                         % interferer center
% legend('target','interferer','interf center')
% 
% % -> fft of interferer
% fft_ran_i = fft(m_time_int, num_samples, 2)/num_samples;       % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
% fft_ran_i = 10*log10(abs( fft_ran_i ).^2);                     % in dB
% 
% % -> fft of interferer + target
% fft_ran_i_t = fft(m_time_int + m_time_targ, num_samples, 2)/num_samples;  % fft_ran( , ,2) along the rows of time_data: fft across the samples in 1°, 2°,... chirps.. Gives fr
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

end