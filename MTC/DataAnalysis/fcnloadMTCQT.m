function [input, D, es] = fcnloadMTCQT(pathname,input,flags,handles,plotflag,ei)
if exist('ei','var') && ~isempty(ei) || flags.status.MC;  MCflag=true; else MCflag=false; end
if MCflag; plotflag = 0; end;  fprintf('Loading MTC data... '); tic

%GET COLLECT DIRECTORY
if ~isfield(input,'MTC') || ~isfield(input.MTC,'rundir') || (exist('rundir','var') && ~isempty(pathname) && ~strcmp(pathname,input.MTC.rundir))
    if exist(pathname,'file')==2
        filename=''; nf=1;
    elseif ~exist('rundir','var') || isempty(pathname)
        [filename,pathname]=uigetfile('/Users/glennjocher/downloads/*.*','Select MTC file:',fcnlastfileloaded,'Multiselect','on'); if pathname==0; fprintf('No file selected ... Done.\n'); return; end
        nf=1; 
        if iscell(filename); nf=numel(filename); end
    end

    if nf==1 && any(strfind(filename,'.mat')) %matfile
        mfname=[pathname filename];
        fprintf('Loading %s... ',mfname);  tic;  load(mfname);  fprintf('Done (%.1fs).\n',toc);
    elseif nf==1 %textfile
        A=fcnloadtextfile([pathname filename],10E6,3);
        if ~isfield(A,'E'); A.E=A.x(:,2);  A.unixTime=A.x(:,1);  A.events=fcnunique(A.E);  A.x=int16(A.x); end
%         if A.loadTime>1
%             mfname = [pathname filename '.mat'];
%             fprintf('Saving ''.mat'' file... '); tic; b=whos('A'); b=b.bytes;  if b>2E9; save(mfname,'A','-v7.3'); else save(mfname,'A','-v6'); end; fprintf('Done (%.1fs).\n',toc); 
%         end %v7.3 required for files > 2Gb
    else %multiple files
        for fi=1:nf
            if any(strfind(filename{fi},'.mat'))
                fprintf('Loading %s... ',[pathname filename{fi}]);  tic;  load([pathname filename{fi}]);  fprintf('Done (%.1fs).\n',toc);
            else
                A=fcnloadtextfile([pathname filename{fi}],10E6,3);
                if ~isfield(A,'E'); A.E=A.x(:,2)+fi*1E6;  A.unixTime=A.x(:,1);  A.events=fcnunique(A.E);  A.x=int16(A.x); end
%                 if A.loadTime>1
%                     mfname = [pathname filename{fi} '.mat'];
%                     fprintf('Saving ''.mat'' file... '); tic; b=whos('A'); b=b.bytes;  if b>2E9; save(mfname,'A','-v7.3'); else save(mfname,'A','-v6'); end; fprintf('Done (%.1fs).\n',toc);
%                 end %v7.3 required for files > 2Gb
            end
            A.x(:,1)=fi;
            if fi==1; B=A; else B.x=[B.x; A.x]; B.rows=B.rows+A.rows; B.E=[B.E; A.E]; B.unixTime=[B.unixTime; A.unixTime]; end;% clear A
        end
        A=B; clear B;
        A.events=fcnunique(A.E);  A.x=int16(A.x);
    end
    if isfield(input.MTC,'A')
        A.pedestals=input.MTC.A.pedestals;
        A.pedestalOutliers=input.MTC.A.pedestalsOutliers;
        A.pruneWindows=input.MTC.A.pruneWindows;
    else
        if eval(A.filename(5:8))==99
            load('pedestals.exp_0099_run_0070.glenn.mat'); %5VOFF PEDESTAL
        else
            %load('pedestals.exp_0001_run_0882.glenn.mat'); %2016 NIST (HV5ON)
            load('pedestals.exp_0002_run_2668.glenn.mat'); %2017 Hawaii HVOFF
            %load('pedestals.exp_0002_run_2689.glenn.mat'); %2017 Hawaii HVON
        end
        %Run 882 - HV on, 5V on
        %Run 883 - HV off, 5V on
        %Run 884 - HV off, 5V off
        A.pedestals=muvb;
        A.pedestalOutliers=pedestalOutliers;
        A.pruneWindows = [];
    end %loads X 71x98304
    A.SRCCHi = MTCpixelID2SRCCH(1:1536);
    A.uhj=load('MTCmapping.mat'); A.uhj=A.uhj.X;
    
    %SCROD AND ASIC OFFSETS
    load MTCoffsets.mat; A.offsets=map; %pixel offsets (1-1536);
    
    %NEW GLENN FORMAT
    input.MTC.SRCCHi = A.SRCCHi;
    input.MTC.rundir = pathname;
    input.MTC.runname = A.pf2;
    [A, C] = fcnprocess1run(A,handles,plotflag);
    input.MTC.A = A;
    input.MTC.C = C;
    input.MTC.events = 1:size(input.MTC.C,1);
    input.eventNumber=1;  handles.GUI.edit1.String=1;
    if plotflag; sca(handles.GUI.axes1); end
end
D = zeros(input.cube.dsp.samples,1536);
E = input.MTC.events;
ei = floorandceil(input.eventNumber+1,1,numel(E)); %event index
if MCflag; input.eventNumber=ei; end;  event=E(ei)-1;
if ~MCflag; handles.GUI.text8.String = sprintf('Event (%g-%g), MTC event %g',[minmax3(E) event]); end


x=input.MTC.C{ei,1}; %event data
es=full(input.MTC.C{ei,2}); %event stats

%nw=min(i+511,input.MTC.A.n*64);
nw = size(x,1);
i=find(any(x,2),1,'first');  j=i:nw;  D(1:numel(j),:)=x(j,:);

%if plotflag;  eventViewer(input, []);  sca(handles.GUI.axes1);  end
fprintf('Done (%.1fs).\n',toc)
end
