% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function fcnrunTSCL()
close all; clc; clear; evalin('base','close all; clc; clear');
%load M2.mat;  input.MTC.A.x=[]; input.MTC.C=[];
%load M2handles.mat


%FOR MTC GEANT MCs
% if ismac %MACBOOK
%     up='/Users/glennjocher/Downloads/';
%     gf={'GEANT/antineutrino_10umStepSizeLimit_1percentNaturalBoron_noNKillerProcess_run0.dat.txt'...
%         'GEANT/1MeV point source.txt'};
%         %'GEANT/EJ-254 1% Boron/gammas/60Co_gammas_100m_cube_0.01mm_step.dat'};
% else %MTCB
%     up=userpath;
%     gf={'GEANTfiles/1MeV point source.txt',...
%         'GEANT/antineutrino_10umStepSizeLimit_1percentNaturalBoron_noNKillerProcess_run0.dat.txt',...
%         'GEANTfiles/EJ-254 1% Boron/gammas/60Co_gammas_100m_cube_0.01mm_step.dat'};
% end
% flags.status.mtc=0;

%FOR MTC NEUTRON SKYMAPS (REAL DATA)
load M431.mat; input.MTC=[];
p='/Users/glennjocher/Downloads/';
%f={'exp_0002_run_2740.glenn','exp_0002_run_2741.glenn','exp_0002_run_2742.glenn'};
[f,p]=uigetfile([p '*.glenn'],'Multiselect','on');
for i=1:numel(f)
    fcnrunTS(fcnloadMTCQT([p f{i}],input,flags,handles,0), GEANT, handles, flags);
end


% %FOR FIBER64 MCs
% close all; clc; clear
% load 'FIBER64_sansGUI.mat';
% flags.status.mtc=0;
% p='/Users/glennjocher/Downloads/GEANT/fiberFiles/';
%      f={'gamma_run_50k_events_0to3MeV.dat.txt'...
%          'neutron_run_50k_events_0to10MeV.dat.txt'};
% fcnrunTS(input, fcnloadGEANT([p(1:end-1) filesep f{1}]), handles, flags)
