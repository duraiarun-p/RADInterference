clc;clear;close all;
%%
%Parameters for chirp
F1=2e6; %Start
F2=3e6; %Stop
TimeDuration=20e-6;%sec

Fs=F1*10; %Sampling 10 times of F1


HalfTime=0.5*TimeDuration; % Cross over time

t=0:1/Fs:TimeDuration;t=t';
Sig=chirp(t,F1,HalfTime,F2);

%%
figure(1),
subplot(121),
% pspectrum(Sig,Fs,'spectrogram','TimeResolution',0.1,'OverlapPercent',99,'Leakage',0.85)
plot(t,Sig);
subplot(122),
pspectrum(Sig,Fs,'spectrogram');