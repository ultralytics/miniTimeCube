function [] = plotAnalysis()
clc; clear; close all 
%MTC RUNS: https://docs.google.com/spreadsheets/d/1LEZYzBG8Ae_jQcYmxaTPxrfccSzXTMRNtyozElfS98k/edit?usp=sharing

%load /Users/glennjocher/Downloads/DATA/MTC/RESULTS7-2016.11.01.mat
load exp_0002_run_2465_subtracted.glenn.mat.fit.mat
%load exp_0001_run_0165.glenn.txt.mat
if ~exist('X','var')
    p=[fcnpathm1(mfilename('fullpath')) 'runAnalysis'];
    a=what(p);  a=a.mat;  n=numel(a);  X=cell(n,1); F=X; P=X; B=X;
    for i=1:n
        if regexp(a{i},'.glenn.')
            b=load([p filesep a{i}]);
            X{i}=b.X;  F{i}=b.F;  P{i}=b.P;  B{i}=b.B;
        end
    end
    i=~cellfun(@isempty,F);  F=F(i);  P=P(i);  X=X(i);  B=B(i);  X=cat(1,X{:}); %#ok<NASGU>
    mfp = [fcnpathm1(mfilename('fullpath')) 'runAnalysis' filesep];
    if numel(P)>30
        s=['RESULTS-' datestr(now,'yyyy.mm.dd')];
        save([mfp s '.mat'],'X','F','P','B');
        gfit(X,[mfp filesep s])
        return;
    end
end
X=cell2mat(X(:,[3 1 2]));
try
    [~,category] = getRunCategories;  %h=fig; histogram(a,'binmethod','integers'); axis tight; h.XTick=0:14;
    c=category(X(:,1));  %plotCategories(c);
catch
    c = ones(size(X,1),1);
end

i=any(c==[1 2 4 5],2);  fprintf('%g category (1,2,4,5) events over %g runs for %.4g days of reactor ON live-time\n',sum(i),numel(unique(X(i,1))),lifetime(X(i,:)));
i=any(c==6,2);  fprintf('%g category 6 events over %g runs for %.4g days of reactor OFF live-time\n',sum(i),numel(unique(X(i,1))),lifetime(X(i,:)));
i=any(c==[1 2 4 5],2);  
X=X(i,:);  ne=sum(i);

%Pe=[-0.41457       2.0597      -3.6936        2.331      0.76008      0.046772]; %Co-60 1 MeV
%Pe=[6.9806e-05   -0.0025071    0.0033318     0.087172       1.2219     0.018164]; %0-10 MeV Kurtis GEANT
%Pe=[-0.00073035     0.011997    -0.076975      0.21981       1.2163   -0.0048373]; %0-10 MeV Point Source
load MTCEcalPoint.mat
X(:,10+3)=F(X(:,10+3));  %635 Pixel MTC ENERGY CALIBRATION
X(:,5+3)=F(X(:,5+3));  %635 Pixel MTC ENERGY CALIBRATION

%CANDIDATES
% %Kurtis PMT Flashing skim cuts 19.6% of events in run 817
% ps=~all(X(:,(1:5)+3)==0 | isnan(X(:,(1:5)+3)),2); %prompt
% ds=~all(X(:,(1:5)+8)==0 | isnan(X(:,(1:5)+8)),2); %delayed
% gt=X(:,3)~=0; %good time
% ptb=X(:,19)>50; %prompt trigger bits
% pe=X(:,5+3)>.3; %prompt energy
% pcw=X(:,10+3+1)>120; %prompt candidate waveforms
% pts=X(:,23)<20; %prompt time sigma
% dts=X(:,35)<20; %delayed time sigma
% dsw=X(:,33)>.1; %delayed sine waves
% %cc=[ps ~ds gt ptb pe pcw pts dts*0+1 dsw*0+1]; %source run
% cc=[ps ds gt ptb pe pcw pts dts dsw]; %neutrino run
% i=all(cc,2); ne=sum(i);  X=X(i,:);
% vpa(mean([cc i])',3)

%c=category(X(:,1)); plotCategories(c); [cu,~,~,nu]=fcnunique(c); d=zeros(15,1); d(cu+1)=nu; 

% %PLOT TIMES
% try
%     load('NBSR.data.7.1.to.10.12.2016.mat')
%     t=unix2now(X(:,3)/1E6);  day=t;
%     ha=fig(3,1,1,2);  x=A(:,1);  s={'NC-3 LINEAR POWER (%)','NC-4 LINEAR POWER (%)','THERMAL POWER (MW)'}; %#ok<NODEF>
%     sca; plot(x,A(:,2),'Display',s{1}); plot(x,A(:,3),'Display',s{2}); datetick; xyzlabel('','POWER (%)'); legend show
%     sca; plot(x,A(:,4),'Display',s{3}); datetick; xyzlabel('',s{3}); legend show;
%     sca; [y,x]=histcounts(day,min(day):1/1440:max(day));  plot(x(2:end),y,'Display','MTC');  datetick('x','mmm'); xyzlabel('','#/minute'); legend show; ha(3).YLim(2)=300;
%     fcntight('xjoint');
% end

%PLOT POINT FITS
load M2.mat; input=eval('input');
load MTCoffsets.mat;  [~,good] = fcnpruninglist;
xp=X(:,(1:5)+3);  xd=X(:,(1:5)+8);
for i=1
    switch i
        case 1; x=xp; s='NEUTRON ANALYSIS 1';
        case 2; x=xd; s='DELAYED ANALYSIS 7';
    end
    x=x(~all(x==0 | isnan(x),2),:);  if ~any(x); continue; end
    ha=fig(1,3,1.1);  
    %ha=fig(2,2,1,1.1);  
    %ha=fig(1,4,1.1,1);
    c=fcndefaultcolors(i);

    nb=100;
    sca; histogram(x(:,4),nb,'facecolor',c); xyzlabel('T (ns)');  title(sprintf('%g %s Events',size(x,1),s)); ha(3).Title.FontSize=16;
    sca; a=x(:,5); histogram(a,linspace(0,max(fcnsigmarejection(a,4,3)),nb),'facecolor',c,'DisplayName','Measured'); xyzlabel('E (MeV)');
    addSimulatedSpectrum(i); %title(sprintf('%g %s Events',ne,s)); ha(3).Title.FontSize=16;
    sca; plot3(x(:,1),x(:,2),x(:,3),'.','color',c); fcnPlotDetector(input,ones(1536,1).*good*.1);  ha(3).CameraViewAngle=7; axis off

%     sca;  nb=400;  zv=zeros(nb,nb);  v=linspace(-66,66,nb);
%     fcnPlotDetector(input,ones(1536,1).*good*.1);  axis off; cla
%     box=input.cube.box;  handles.detectoroutline = plot3(box.x(:), box.y(:), box.z(:),'Color',[.7 .7 .7]);
%     [N,xe,ye] = histcounts2(x(:,1),x(:,2),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,ye,zv-67,N,'edgecolor','none');
%     [N,xe,ye] = histcounts2(x(:,2),x(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(zv+67,xe,ye,N,'edgecolor','none');
%     [N,xe,ye] = histcounts2(x(:,1),x(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,zv-67,ye,N,'edgecolor','none');  ha(4).CameraViewAngle=7;
        
    fcnmarkersize(1.5); fcntight; fcntight('csigma')
end


% %PLOT ANGLE DISTRIBUTION
% fig(1,3,1.5);
% for i=1:3
%     switch i
%         case 1; x=xp; s='PROMPT';
%         case 2; x=xd; s='DELAYED';
%         case 3; x=xd-xp; s='DELAYED - PROMPT';
%     end
%     if ~any(x); break; end
%     sca; d=fcnvec2uvec(x(:,1:3));  ea=fcnelaz(d);  m=2;
%     c = histcounts2(ea(:,1)*r2d,ea(:,2)*r2d,linspace(-90,90,180*m+1),linspace(-180,180,360*m+1));  %c=x(2:end)/2+x(1:end-1)/2;
%     [xx,yy] = ndgrid(linspace(-90,90,180*m),linspace(-180,180,360*m));
%     [rm,cm] = fcnmapprojection(xx,yy,'winkeltripel');  gain=min(cm(:,1))./cm(:,1);  cg=c.*gain;
%     pcolor(cm,rm,cg); shading flat;
%     fcnplotdetectorprojection(input,ones(1536,1),ones(1536,1)*200,[0 0 0],'winkeltripel','');  title(sprintf('%g %s',ne,s));
%     h=gca; h.CameraViewAngle=7.5; h.Position(1)=h.Position(1)-.02;
% end
% fcnfontsize(12); fcntight('csigma')

%PLOT STATISTICS
x=X(:,14:end);     s={'Candidate Waveforms','Amplitude','Width','Time','Integral','Triggers'};  t={' \mu',' \sigma'};
fig(2,6);
for j=1:2
    for i=1:6
        xi=fcnsigmarejection(x(:,i+(j-1)*6),6,2);
        sca; histogram(xi,75);  xlabel([s{i} t{j}])
        if j==1 && i==3; title(sprintf('%g PROMPT Events',ne)); ha(3).Title.FontSize=16; end
    end
end
fig(2,6);
for j=1:2
    for i=1:6
        xi=fcnsigmarejection(x(:,i+(j-1)*6+12),6,2);
        sca; histogram(xi,75,'FaceColor',fcndefaultcolors(2));  xlabel([s{i} t{j}])
        if j==1 && i==3; title(sprintf('%g DELAYED Events',ne)); ha(3).Title.FontSize=16; end
    end
end

%DT
fig(1,3,1.1,1);  v=linspace(0,10,100);
sca; histogram(X(:,7)/1E3,v); xyzlabel('T_{PROMPT} (\mus)');
sca; histogram(X(:,12)/1E3,v,'FaceColor',fcndefaultcolors(2)); xyzlabel('T_{DELAYED} (\mus)','','',sprintf('%g ANALYSIS 7 Events',ne));
sca; histogram(X(:,12)/1E3-X(:,7)/1E3,v,'FaceColor',fcndefaultcolors(3)); xyzlabel('DT (\mus)');

%WAVEFORMS VS ENERGY
%hist211(X(:,21),xd(:,5),{linspace(0,200,200),linspace(0,1,200)});

%PROMPT VS DELAYED ENERGY
hist211(X(:,8),X(:,13),{linspace(0,6,100),linspace(0,.5,100)}); xyzlabel('E_{PROMPT} (MeV)','E_{DELAYED} (MeV)'); fcnmarkersize(0.1);

end


function []=addSimulatedSpectrum(i)
load simulatedNuebarSpectrum.mat
X=NS{i}{1};  Y=NS{i}{2}; %#ok<USENS>

h=gca; x=linspace(h.XLim(1),h.XLim(2),1000);  y=interp1(X,Y,x,'linear');
y=y./max(y)*max(h.Children.Values(:));
plot(x,y,'k','Display','Simulated','LineWidth',2); legend show
end


function gfit(X,filename)
T=cell2mat(X(:,[3 1 2]));  T(:,3)=T(:,3)/1E6;
s={'runNumber'
'eventNumber'
'unixTime'
'X'
'Y'
'Z'
'T'
'E'
'X'
'Y'
'Z'
'T'
'E'
'sumTriggers'
'meanAmplitude'
'meanWidth'
'meanTime'
'meanIntegral'
'sumWaveforms'
'sigmaTriggers'
'sigmaAmplitude'
'sigmaWidth'
'sigmaTime'
'sigmaIntegral'
'sigmaWaveforms'
'sumTriggers'
'meanAmplitude'
'meanWidth'
'meanTime'
'meanIntegral'
'sumWaveforms'
'sigmaTriggers'
'sigmaAmplitude'
'sigmaWidth'
'sigmaTime'
'sigmaIntegral'
'sigmaWaveforms'};

fid = fopen([filename '.gfit'],'w');
fprintf(fid,'# Ultralytics LLC\n# MTC fit results by glenn.jocher@ultralytics.com on %s\n',datestr(now,'yyyy.mm.dd-HH.MM.SS'));
fprintf(fid,'# For more information see: https://docs.google.com/spreadsheets/d/1SWokfpsM5a0YlEPJ6AtCLcp-xwXBPmOiasHeNz7B4uw/edit#gid=700280350\n#\n# ');
fprintf(fid,['%10s%15s%30s' repmat('%10s',1,34) '\n'],char(s)');
fprintf(fid,['%10g%15g%30.20g' repmat('%10g',1,34) '\n'],T');
fclose(fid);
end

function days=lifetime(X)
runs=unique(X(:,1));  n=numel(runs);  T=X(:,3)/1E6;  i=T~=0 & ~isnan(T); T=T(i); X=X(i,:);  t=0;
for i=1:n
    ti=T(X(:,1)==runs(i));
    if any(ti); t = t + (max(ti)-min(ti)); end
end
days=t/86400;
end

function []=plotCategories(c)
h=fig(1,1,'15cm');
d=c(any(c==[1 2 4 5],2));
histogram(c,'binmethod','integers','FaceColor',[.7 .7 .7],'DisplayName','All');
histogram(d,'binmethod','integers','FaceColor','g','DisplayName','\nu Candidate Categories');
 axis tight; h.XTick=0:14; xyzlabel('MTC Run Category','Events','',sprintf('Analaysis 7 - %g Green Category Events',numel(d))); legend show
end