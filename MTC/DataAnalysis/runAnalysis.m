% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = runAnalysis()
close all; clc; clear all %#ok<CLALL>
load M2.mat;  input.MTC.A.x=[]; input.MTC.C=[]; %#ok<*STRNU>
%pathname = uigetdir('/Users/glennjocher/downloads/','Select Directory:');
%pathname = '/Users/glennjocher/Downloads/MTC';

pathname=cell(10,1);
for i=1:10
    pathname{i} = uigetdir('/Users/glennjocher/downloads/','Select Directory:');  if pathname{i}==0; break; end
end; pathname=pathname(1:i-1);  %save gP.mat pathname
%pathname = {'/home/kurtisn/skimmedData','/scratch/skimmedData','/data/skimmedData','/data2/skimmedData','/data3/skimmedData'};

gP={};
for i=1:numel(pathname)
    gP=cat(1,gP,genpathcell(pathname{i}));
end
np=numel(gP);
mfp = [mfilename('fullpath') filesep];  if ~exist(mfp,'dir'); mkdir(mfp); end;  input.mfp=mfp;


F=cell(1E5,1); P=F; B=F; startclock=clock;
for i=1:np
    D=dir(gP{i});
    for j=1:numel(D)
        if isempty(D(j).bytes); D(j).bytes=0; end
    end
    %     fi=find(~cellfun(@isempty, strfind({D.name},'.glenn')) & ...
    %         cellfun(@isempty, strfind({D.name},'.mat')) & ...
    %         [D.bytes]>1E6);  nf=numel(fi);
    fi=find(~cellfun(@isempty, strfind({D.name},'.glenn.mat')) & ...
        [D.bytes]>1E6);  nf=numel(fi);
    if nf>0
        j=ones(nf,1)*i;
        F{i}={D(fi).name}';
        P{i}={gP{j}}'; %#ok<CCAT1>
        B{i}={D(fi).bytes}';
    end
end
F=cat(1,F{:});  P=cat(1,P{:});  B=cat(1,B{:});
[~,i]=unique(F);  F=F(i);  P=P(i);  B=B(i);  n=numel(F);

% %GREEN RUNS ONLY
% j=false(n,1);  [~,category] = getRunCategories;
% for i=1:n
%     a=regexpi(F{i},'_run_(.{1,6}).glenn','tokens');  run=str2double(a{1}{1});
%     goodrun=any(category(run)==[1 2 4 5 6]);
%     j(i)=~exist([mfp F{i} '.mat'],'file') && goodrun;
% end; F=F(j); P=P(j); B=B(j);

[~,i]=sort(cell2mat(B));  F=F(i);  P=P(i);  B=B(i); %order by size

n=numel(F);  notdone=true;
%while notdone
    %try
        %delete(gcp('nocreate'));  parpool(4,'SpmdEnabled',false);
        for i=1:n %parfor
            matname=[mfp F{i} '.fit.mat'];
            if ~exist(matname,'file')
                try
                    X=oneFile(input,[P{i} filesep F{i}]);
                    parsave(matname,X,F{i},P{i},B{i});
                catch ME
                    parsave([mfp F{i} '.fit.fail'],ME,F{i},P{i},B{i});
                end
            end
        end
    %catch
    %    notdone=true;
    %    continue
    %end
    %notdone=false;
%end
e = etime(clock,startclock);
fprintf('Finished %g files in %g paths in %.1fhrs (%.0fs).\n\n',n,np,e/3600,e)

plotAnalysis;
end


function [C, Arows, Aftell]=oneSegment(input,filename,maxrows,skiprows)
A=fcnloadtextfile(filename,maxrows,skiprows);  Arows=A.rows;  Aftell=A.ftell;

A.unixTime=A.double;  A.double=[];  A.x=int16(A.x);  %A.E=A.x(:,2);
A.pedestals=input.MTC.A.pedestals;  A.pedestalOutliers=input.MTC.A.pedestalOutliers;  A.SRCCHi=input.MTC.A.SRCCHi;  A.uhj=input.MTC.A.uhj;  A.offsets=input.MTC.A.offsets;
input.MTC.A=[]; input.MTC.A.offsets=A.offsets;  %free up space
pid = int16(MTCSRCCH2pixelID(A.x(:,3:6)));  nw=512;  RW=A.x(:,7);  W=A.x(:,8);  i=W<RW;  Wa=W-RW+1;  Wa(i)=Wa(i)+nw;  clear RW

[~, i] = fcnpruninglist(pid); 
%i = i & Wa<426 & Wa> 418; %FOR NEUTRON RUNS ONLY!!

pid = pid(i);  if isempty(pid); C=[]; return; end
W=W(i);  V=A.x(i,10:73)';  Wa=Wa(i);  A.E=A.E(i);  A.unixTime=A.unixTime(i);  A.x=A.x(i,:);  % Window, nx64 Voltages, Scrod, WindowAdjusted, Event
V = subtractpedestals(A.pedestals,A.pedestalOutliers,pid,W,V);  A.pedestals=[]; A.pedestalOutliers=[];  clear W

i = ~all(V==1600);  pid=pid(i);  V=V(:,i);  Wa=Wa(i);  A.E=A.E(i);  A.unixTime=A.unixTime(i);  A.triggers=A.x(i,9);  A.x=[];
V=removeSpikes(V);
A.nw=double(max(Wa)-min(Wa)+1); %window adjusted, fwi=first window index
A.rows=numel(A.E);  [A.events, A.ei] = fcnunique(A.E);  A.ei=[A.ei; A.rows+1];  ne=numel(A.events);  clear E %event start indices
A.unixTime=A.unixTime(A.ei(1:end-1));
V=V';

i=regexpi(A.filename,'run_');  runNumber=str2double(A.filename(i+4:i+7));

ne=min(100E3,ne);  A.fwi=zeros(ne,1);  Ca=cell(ne,1); Cb=Ca; Cc=Ca;  tic
waitbarParfor(ne);
delete(gcp('nocreate'));  parpool(4,'SpmdEnabled',false);
parfor i=1:ne %parfor
    fprintf('%s\n',waitbarParfor);
    %fprintf('\n%.0f/%.0f  ',i,ne)
    
    [C1, D1]=oneEvent(A,Wa,pid,V,i,'prompt');
    [C2, D2]=oneEvent(A,Wa,pid,V,i,'delayed');
    op=fitsPrepare(C1,C2);  op(1).D=D1;  op(2).D=D2;
    results = fcnantineutrino(input,op,[],[],0);
    
    Ca{i}=results.xhat;
    Cb{i}=[eventStats(C1) eventStats(C2)];
    Cc{i}=[runNumber double(A.events(i)) double(A.unixTime(i))];
end; fprintf('\n Done. (%.1fs)\n',toc);
C=cell(ne,3);  C(:,1)=Ca;  C(:,2)=Cb;  C(:,3)=Cc;
end


function C=oneFile(input,filename)
maxrows=100E3;
[C, Arows, Aftell]=oneSegment(input,filename,maxrows,3);
if Arows==maxrows %need more segments
    for i=1:100
        [C1, Arows, Aftell]=oneSegment(input,filename,maxrows,Aftell+1);
        C=cat(1,C,C1);
        if Arows<maxrows; break; end
    end
end
end


function output=fitsPrepare(C1,C2)
output(1)=fitPrepare(C1);
output(2)=fitPrepare(C2);
    function A=fitPrepare(C)
        D=zeros(256,1536);
        j=find(C(:,1));
        A = oneBookend([j C(j,2:end)],C,D,[],0,1);
    end
end


function [mus,mu,s] = eventStats(a)
k=find(a(:,4)~=0);
[~,j]=fcnsigmarejection(a(k,4),2,2);  a(k(~j),4)=0; %reject bad times

a(:,1)=sum(a(:,4)~=0);
a(:,6)=sum(a(:,6)~=0);
b=a(a(:,4)~=0,:);
mu=nanmean(b,1);
s=nanstd(b,[],1);
mus=[mu s];
end


function [C, D]=oneEvent(A,Wa,pid,V,i,pdstr)
j = A.ei(i):A.ei(i+1)-1;
C=zeros(1536,6); D=[];

switch pdstr
    case 'prompt'
        if (max(Wa(j))-min(Wa(j)))<3; return; end
        j=firstfour(Wa(j),j); %FIRST 4
    case 'delayed'
        if (max(Wa(j))-min(Wa(j)))<7; return; end
        j=lastfour(Wa(j),j); %LAST 4
end

Wj=Wa(j);  A.fwi(i)=min(Wj);  nw=max(Wj)-min(Wj)+1;  nj=numel(j);

xa = zeros(nw*64,1536,'int16');
k = sub2ind([64 nw 1536],ones(nj,1),Wj-min(Wj)+1,pid(j));
xa(k+(0:63)) = V(j,:); %fig; plot(xa)
tv=accumarray(pid(j),A.triggers(j),[1536,1]); %trigger vector

[j, a, w, t, integral, D] = fcnpulsewidth(single(xa),[.45 .50 .55],[50 90000],[8 64],1600,'fraction');  k=t~=0;

%APPLY CALIBRATIONS
a(j) = a(j)./A.offsets(j,2);                            a(isnan(a) | a<0)=0;
t(k) = t(k) - A.offsets(k,4) + A.fwi(i)*64*.360;        t(isnan(t))=0;
integral(j) = integral(j)./A.offsets(j,5);              integral(isnan(integral))=0;

%OUTPUT
C = [k, a, w, t, integral, tv]; %[candidatesJocher, amplitude, width, time, integral, triggersKurtis]
end


function B=genpathcell(pathname)
A = genpath(pathname);
ia=1; fA=find(A==':'); na=numel(fA);  B=cell(na,1);
for i=1:na
    ib=fA(i);
    B{i}=A(ia:ib-1);
    ia=ib+1;
end
end


function parsave(fname,X,F,P,B) %#ok<INUSD>
save(fname,'X','F','P','B')
end