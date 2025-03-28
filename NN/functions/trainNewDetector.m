% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function []=trainNewDetector()
evalin('base','clc; clear all; close all;'); startclock=clock;

%0-SAVE DETECTOR CONFIGURATION
%handles.GUI=[];  handles.GUI.realdataflag=false;  handles.GUI.FiberFlag=true;  save FIBER64_1L_sansGUI.mat


%1-GATHER TRAINING DATA
load 'FIBER64_1L_sansGUI.mat'; input.a=[];
fname=fcnrunTS(input,GEANT,handles,flags);


%2-TRAIN NETWORK AND SAVE
trainFiber64({fname})
return

%3-PARTICLE TRADE STUDY IN DETECTOR
input.NN.fiberfit = load(['NN.' fname]);
%input.NN.fiberfit = load('NN.TS.01 FiberPointSource - 3L SENSL J60035 5V 128ch BCF10 50k.mat');
p='/Users/glennjocher/Downloads/DATA/GEANT/fiberFiles/';
f={'gamma_run_50k_events_0to3MeV.dat.txt'...
         'neutron_run_50k_events_0to10MeV.dat.txt'};

fcnrunTS(input, fcnloadGEANT([p f{2}]), handles, flags, 1E6)

e = etime(clock,startclock);  fprintf('Finished in %.1fhrs (%.0fs).\n\n',e/3600,e)