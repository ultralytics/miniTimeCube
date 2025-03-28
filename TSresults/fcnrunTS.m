% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function filename=fcnrunTS(input,GEANT,handles,flags,nj)
try closeallexcept(handles.GUI.figure1);  flags.status.mtc=ischecked(handles.GUI.realdataflag); end %#ok<TRYNC>

%USER INPUTS --------------------------------------------------------------
if ~exist('nj','var'); nj=100000; end;  if flags.status.mtc;  nj = numel(input.MTC.events); end %number of MCs to run
vnj = repmat(GEANT.ve(:),[ceil(nj/GEANT.ne), 1]);  vnj = vnj(1:nj);
pidname = GEANT.particleName;    %'Antineutrino','Muon','Gamma','Neutron','FiberTrain','FiberNeutron','FiberGamma'
if flags.status.fibers;  pidname=['Fiber' pidname];  end
if flags.status.mtc; pidname='MTCEvents';  vnj=0:(nj-1);  end
%pidname = 'Neutron';

%TRADE STUDY VECTOR -------------------------------------------------------
tsv = 1;
%[GEANTfiles, GEANTpath] = uigetfile('/Users/glennjocher/downloads/*.txt','Select GEANT files:','Multiselect','on');
% GEANTfiles={'gamma_run_50k_events_0to3MeV.dat.txt','neutron_run_50k_events_0to10MeV.dat.txt'};  GEANTpath='/Users/glennjocher/Downloads/DATA/GEANT/fiberFiles/';
% if exist('GEANTfiles','var')
%     tsv=1:numel(GEANTfiles);
%     if strcmp(pidname,'FiberPointSource'); disp('WARNING, NEED TO LOAD EITHER NEUTRON OR GAMMA WITH FIBERS!!'); return; end
% end                                                       
%mex('/Users/glennjocher/Google Drive/MATLAB/neutrinos/nViewGUI/C/photontransportcV4.cpp')

%CREATE NAME FOR MAT FILE
[njp,nja]=numberPrefix(nj);
filename = sprintf('TS.25ns %s - %s %s %.0fch %s %.0f%s.mat',pidname,input.cube.prettyvolume,input.cube.prettyname,input.cube.pixels, input.Material(1).name, nja,njp);
%filename = sprintf('MTC.%s.v3 %s - %.0fch %.0f%s.mat',input.MTC.A.filename(14:17),pidname,numel(fcnpruninglist),nja,njp);
%[filename,pathname] = uiputfile('*.mat','Create a new file to save MC results to:',[cd '/TSresults/' filename '.mat']);  if filename(1)==0; return; end
if ~exist('pathname','var');  pathname=fcnpathm1(mfilename('fullpath')); end
fprintf('Starting ''%s''...\n',filename)

%PREALLOCATE --------------------------------------------------------------
timeOffset = rand(nj,1)*25;
energyGain = ones(nj,1);
flags.status.vertexonwall=1;
switch pidname
    case {'FiberNeutron','FiberGamma'}
        nx = 29;  nt = 29;
    case {'Antineutrino','Neutrino'}
        nx = 10;  nt = 10;  flags.status.vertexonwall=0;
    case 'FiberMuon'
        nx = 5;  nt = 5;
    case 'Muon'
        nx = 4;  nt = 4;
    case 'FiberPointSource'
        if flags.status.lappd
            nx = 1024;  nt = 8;  %LAPPD/NTC: nx=1024 nt=8
        else
            nx = 512;   nt = 4;  %FIBER1: nx=512 nt=4;
        end
        %energyGain = min(random('unif',0,1,[nj 1]),30);          %energyGain(:)=0.2; %diffuse tests
        energyGain = min(random('exp',0.25,[nj 1]),30/4);
        flags.status.vertexonwall=0;
    case 'Gamma'
        nx=5;  nt=5;
    case 'Neutron'
        nx=28; nt=28;
    case 'MTCEvents'
        nx=5;  nt=5;
    case 'PointSource'
        nx=5;  nt=5;  flags.status.vertexonwall=0;
        energyGain = rand(nj,1)*1;
    case 'gammaStats'
        nx=7;  nt=7;
    otherwise
        nx=5;  nt=5;
end

ni                      = length(tsv);
xhat                    = zeros(nj,nx,ni,'single'); %estimates [xyzet]
xtrue                   = zeros(nj,nt,ni,'single'); %truth
vertex                  = zeros(nj,10,ni,'single'); %[p0 t0 e0 v0 energyGain triggers]
collectedPhotonCount    = zeros(nj,2, ni,'single'); %'uint16');

%OPTIMIZER OPTIONS --------------------------------------------------------
input.optimizer.psoptions2 = psoptimset('TimeLimit',3600,'Vectorized','on','CompletePoll','on','MeshAccelerator','on','Display','off','InitialMeshSize',1,'PenaltyFactor',100,'PollMethod','GPSPositiveBasis2N');
input.MTC.C(:,1,:)={[]};

%RUN MONTE CARLOS ---------------------------------------------------------
flags.status.MC=1;  rng(0);  seeds=round(rand(nj,ni,'single')*1E9);  startclock=clock; if ~exist('tsv','var'); tsv=1; end
for i = 1:ni
    tsvi=tsv(i);  input.tsv=tsvi;  startTS=clock;
    %flags.update.detectorgeometry=1;  input=fcnDefinePixels(input,flags,handles);  flags.update.detectorgeometry=0;
    if exist('GEANTfiles','var') && ~strcmp(GEANTfiles{i},GEANT.filename);  GEANT=fcnloadGEANT([GEANTpath GEANTfiles{i}]);  end
    MC.particleName{i}=GEANT.particleName;
        
    waitbarParfor(nj);
    delete(gcp('nocreate'));  parpool(4,'SpmdEnabled',false);
    parfor j = 1:nj %parfor here
        fprintf('TS %g/%g: %s\n',i,ni,waitbarParfor);  
        %fprintf('TS %g/%g : %s\n',i,ni,num2str(j));
        rng(seeds(j,i));
        [input1,output,G1,photons,PE] = fcn1Event(input,flags,handles,GEANT,[tsvi timeOffset(j) energyGain(j) vnj(j)]);
        
        %ESTIMATION -------------------------------------------------------
        results=[];
        try
            switch pidname
                case 'Neutron'
                    results = fcnfastneutron(input1,output,handles,PE,G1,0);
                    %results = fcndisambiguate(input,output,handles,photons,flags,0);
                case 'Muon'
                    results = fcnmuon(input1,output,handles,photons,G1,0);
                case {'FiberNeutron','FiberGamma'}
                    %results = fiberNeutronFitter(input1,output,flags,PE,G1,handles,0);
                    results = fcnNeutronExits(input1,output,flags,PE,G1,handles,0);
                case 'FiberMuon'
                    results = fiberNeutronFitter(input1,output,flags,PE,G1,handles,0);
                case 'FiberPointSource'
                    results = gatherTrainingData(input1,output,PE,G1);
                case {'Antineutrino','Neutrino'}
                    results = fcnantineutrino(input1,output,handles,G1,0);
                case {'MTCEvents','Gamma','PointSource'}
                    PE.n = sum(output(1).N); PE.pixel=output(1).pid;
                    results = runMLpoint(input1,output,handles,G1,0);
                case {'gammaStats'}
                    results = gammaStats(G1);
            end
        catch
           %save(sprintf('FAIL.%g.mat',j),'input1','output','handles','PE','G1'); %NO SAVE WITH PARFOR!!
        end
        try
            xhat(j,:,i)  = results.xhat;
            xtrue(j,:,i) = results.true;
        catch
            sprintf('WARNING: fcnrunTS.m FAILURE for event %g',j)
        end
        
        if PE.n>0; triggers = numel(fcnunique(PE.pixel));  else;  triggers = 0; end % number of triggers
        vertex(j,:,i)                 = [G1.p0 timeOffset(j) G1.uake1(1) G1.uvec(1,:) energyGain(j) triggers];
        tsigma = nanstd(fcnsigmarejection(output(1).t(output(1).t~=0)));
        collectedPhotonCount(j,:,i)   = [PE.n tsigma];
    end
    delete(gcp('nocreate')); 
    
    e=etime(clock,startTS);  MC.elapsedHours(i)=e/3600;
    fprintf('%.4g TS point (%.0f Events) Finished in %.1fhrs (%.0fs).\n\n',tsv(i),nj,e/3600,e)
end
rng(0);
clear input1

%ALLOCATE TO MC STRUCTURE -------------------------------------------------
MC.xhat                     = xhat;
MC.xtrue                    = xtrue;
MC.vertex                   = vertex;
MC.collectedPhotonCount     = collectedPhotonCount;
MC.GEANT.files              = GEANT.files;

%PRINT WORKSPACE REPORT ---------------------------------------------------
e = etime(clock,startclock); %elapsed
fprintf('\n\n%.0f MC Runs Completed in %.1fhrs (%.0fs), %.3fs per run.\n',ni*nj,e/3600,e,e/ni/nj)
fprintf('Saving ''%s''\nin folder ''%s''...\n',filename,pathname)

%SAVE RESULTS TO A MAT FILE
input.MTC.A.pedestals=[];  input.MTC.A.pedestalOutliers=[];  input.MTC.A.x=[];  input.MTC.C=[];
MC.FileName = filename;
MC.date = date;
s=whos('MC');
if s.bytes>1.9E9 % MC > 2 Gb
    save([pathname filename],'MC','tsv','input','-v7.3')
else
    save([pathname filename],'MC','tsv','input')
end
delete(gcp('nocreate'));

%PLOTTING -----------------------------------------------------------------
%fcnplotneutronTS(input,MC,tsv)
%fcnplotTSresults(input,MC,tsv)
end


function results = gammaStats(G1)
results.xhat = zeros(1,7);
results.true = zeros(1,7);

x=[rangec(G1.p1,G1.p2) G1.de]; i=find(G1.pid==11); j=find(G1.tid==G1.tid(i(1)));

[~, lv1] = fcntid2ctid(G1.tid(j(1)),G1.tid,G1.ptid);  k=[j; find(lv1)];
results.true = [sum(x(i,:),1) sum(x(j,:),1) sum(x(k,:),1) G1.ke1(1)-G1.ke2(1)];

%PLOT
%fig(4,2);  T=MC.xtrue; for i=1:7; sca; fhistogram(T(:,i),100); end
end
