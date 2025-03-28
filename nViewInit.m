% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

%clc
ha = handles.GUI.axes1;
rng(1);

if ~exist('GEANT','var') %if no MeV table is ed, must be startup. Load MeV Table
    %ADD PATHS ------------------------------------------------------------
    cp = fcnpathm1(mfilename('fullpath')); cd(cp); addpath(genpath(cp))
    cp = fcnpathm1(cp); addpath(cp)
    input.directory = strrep(pwd,'\',filesep);
    
    %SETUP GUI LOADING IMAGE ----------------------------------------------
    ht = text(0,0,0,'Loading GEANT Data...','FontSize',25,'FontName','Arial','HorizontalAlignment','center','Color',[0 160 217]/255);
    if ~exist('output','var') %first time opening
        GEANT = fcnloadGEANT('1MeV point source.txt',1E4);
    else
        GEANT = fcnloadGEANT('',1.2E7); %about 1.2 Gb GEANT .txt file
    end
    G1=[];  photons=[];  PE=[];  output=[];
    delete(ht)
end
hc = get(0,'children');
delete(hc(hc~=handles.GUI.figure1))


%DEFINE SOME INPUTS -------------------------------------------------------
cd(input.directory)
input = fcninittables(input);
input.MTC = [];
%input.NN.fiberfit = load('NN.TS.01 2.2mAL 250keVexp 5mm FiberPointSource - 3L SENSL J60035 5V 128ch BCF10 40k.mat');
%input.NN.fiberpid = load('NN.PATTERN.TS FiberNeutron - 2L SENSL J60035 5V 50pc 512pix EJ-254 1M.mat');

input.neutron = load('dEfit.mat');

%OPTIMIZER OPTIONS --------------------------------------------------------
input.optimizer.gamma = load('gammaKleinNishinaF.mat');
try
    input.optimizer.psoptions1 = psoptimset('Vectorized','on','Display','off');
    input.optimizer.psoptions2 = psoptimset('Vectorized','on','CompletePoll','on','MeshAccelerator','on','Display','off');
    input.optimizer.psoptions3 = psoptimset('Vectorized','on','CompletePoll','on','MeshAccelerator','on','Display','off','PollMethod','GPSPositiveBasis2N');
    input.optimizer.options1 = optimset('Display','off','Algorithm', 'interior-point', 'GradObj','off','Hessian', 'bfgs','SubproblemAlgorithm', 'ldl-factorization');
    input.optimizer.options2 = optimset('Display','off');
catch %#ok<CTCH>
    fprintf('\nWARNING: MATLAB GLOBAL OPTIMIZATION TOOLBOX NOT DETECTED\n')
end

%INITIALIZE UPDATE FLAGS --------------------------------------------------
initFlags

%SET ALL AXES HANDLES -----------------------------------------------------
handles.uvec = [];
handles.event = [];
handles.upid=[];
handles.upidsp=[]; %scintillation photons
handles.upidcp=[]; %cherenkov photons
handles.detector=[];
handles.gridanode=[];
handles.stripanode=[];
handles.segments=[];
handles.detectoroutline=gobjects;
handles.detectorfacelabels=[];

%SETUP PLOT ---------------------------------------------------------------
view(ha,150,20);
xyzlabel(ha,'X (mm)','Y (mm)','Z (mm)'); fcnfontsize(ha,10)
set(ha,'ZDir','reverse','YDir','reverse','CameraViewAngle',6.75);
axis(ha,'equal','vis3d','manual')
colormap(ha,parula(256))

%CLEAR VARIABLES FROM WORKSPACE -------------------------------------------
clear ha ht hc str LoadImage cp

