% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = plotNTCfiber(files)
clc; clear;
if nargin==0
    files=uigetfile('*.mat','MultiSelect', 'on');  if ischar(files); files={files}; end;
    %files={'TS FiberPointSource - 1L SENSL J60035 5V 128ch EJ-254 200k.mat'};
end
nf=numel(files);
for fi=1:nf
    f = files{fi};  fprintf('Loading ''%s''... ',f);  load(f);  fprintf('Done.\n');
    inputs=MC.xhat;  targets=MC.xtrue;  clear input tsv MC
    
    vi{1}=1:256; vi{2}=vi{1}+256;  a1=max(inputs(:,vi{1}),[],2);  a2=max(inputs(:,vi{2}),[],2);  goodAmplitude=a1>.002 & a2>.002;  %goodRatio=max(a1./a2,a2./a1)<2.5;
    i=find(~all(inputs==0,2) & ~all(targets==0,2) & ~any(isnan(inputs),2) & goodAmplitude & targets(:,3)>10);  i=i(1:end); inputs=inputs(i,:);  targets=targets(i,:);  n=size(inputs,1);  tbias=zeros(n,3);  numel(i)
    
    t=zeros(n,2);  a=t; w=t; int=t;  %F=lowPassFilter(100E6,200E6,5E9);
    for i=1:2
        I = inputs(:,vi{i});
        %I = filter(F,I')';  nb=2^12;  I=round(I*nb)/nb; %12-bit digitized
        H = filterV(I(1,:)'); I=filter(H,I')';
        [~,a(:,i),w(:,i),t(:,i),int(:,i)]=fcnpulsewidth(I',.5,[0 1000],[15 50],1E6,'fraction');  inputs(:,vi{i}) = I;
    end;  t(t==0)=inf;
    i=all(w>15,2);% & targets(:,1)>-440 & targets(:,1)<320;
    xa=targets(i,1);  dt = t(:,2)-t(:,1);  A.filename=f;
    
    
    clear h;
    if any(diff(xa)); staticSource=false; else; staticSource=true; end
    if staticSource
        ha=fig(2,4,1); nb=100;
        sca; b=fcnsigmarejection(t(i,1:2)); xa=linspace(minnonzero(b),max3(b),nb)*.2; h(1)=fhistogram(t(i,1)*.2,xa); h(2)=fhistogram(t(i,2)*.2,xa); xyzlabel('(ns)','','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
        sca; b=fcnsigmarejection(a(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(a(i,1),xa); fhistogram(a(i,2),xa); xyzlabel('(V)','','','amplitude')
        sca; xa=linspace(20,60,nb); fhistogram(w(i,1),xa); fhistogram(w(i,2),xa); xyzlabel('(bins)','','','width')
        sca; b=fcnsigmarejection(int(i,1:2),4,2); xa=linspace(0,max3(b),nb); fhistogram(int(i,1),xa); fhistogram(int(i,2),xa); xyzlabel('(V)','','','integral')
        sca; h=fhistogram((t(i,2)-t(i,1))*.2,nb); xyzlabel('dt (ns)','','',''); h.FaceColor=[.7 .7 .7];
        sca; h=fhistogram(a(i,1)./a(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
        sca; h=fhistogram(w(i,1)./w(i,2),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
        sca; h=fhistogram(int(i,1)./(int(i,2)),nb); xyzlabel('ratio','','',''); h.FaceColor=[.7 .7 .7];
        for j=1:8; legend(ha(j),'show'); legend(ha(j),'boxoff'); ha(j).YTick=[]; ha(j).YColor=[1 1 1]; end
        fcntight; text(ha(2),-.04,ha(2).YLim(2)*1.14,f);  fcnlinewidth(1.5);  fcnfontsize(17);
        export_fig('-a1','-q90','-r250',[A.filename '.png'])
    else
        ha=fig(2,4,1);  type='profile';  nb=round(diff(minmax3(xa))/10);  %nb=100;
        YL=[2 12; 0 .03; 15 40; 0 .8; -10 10;0 5; .6 1.6;  0 5];
        sca; for j=1:2; h(j)=hist2(xa,t(i,j)*.2,nb,type); end; xyzlabel('','(ns)','','time'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best');
        sca; for j=1:2; b=a(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,a(i,1)+a(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','amplitude'); legend show; legend boxoff; legend('Location','northwest'); delete(h)
        sca; hist2(xa,w(i,1),nb,type); hist2(xa,w(i,2),nb,type); xyzlabel('','(bins)','','width')
        sca; for j=1:2; b=int(i,j); h=hist2(xa,b,nb,type); exp1fit(h,xa,b); end;  h=hist2(xa,a(i,1)+a(i,2),nb,type); h.DisplayName='sum'; xyzlabel('','(V)','','integral'); legend show; legend boxoff; legend('Location','northwest'); delete(h)
        sca; b=dt(i,:)*.2; h=hist2(xa,b,nb,type); sx=''; sy=''; try h.Color=[.7 .7 .7]; sy=h.DisplayName; h=hist2(b,xa,nb,'profile'); sx=h.DisplayName; delete(h); end;  F=fit(xa,b,'poly1'); hl=plot(xa,F(xa),'Display',sprintf('%.1f mm/ns',1./abs(F.p1))); xyzlabel(sprintf('X (mm) %s',sx),sprintf('(ns) %s',sy),'','dt'); legend(hl,'Location','northwest'); xsigma=sx; tsigma=sy; slope=1./abs(F.p1);
        sca; b=a(i,1)./a(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
        sca; h=hist2(xa,w(i,1)./w(i,2),nb,type); xyzlabel('X (mm)','','','ratio'); try h.Color=[.7 .7 .7]; end
        sca; b=int(i,1)./int(i,2); h=hist2(xa,b,nb,type); F=fit(xa,b,'exp1');  h.DisplayName=sprintf('\\lambda = %.0f mm',1/abs(F.b)); legend show; legend('Location','northwest'); sx=''; try h.Color=[.7 .7 .7]; h=hist2(b,xa,nb,type); sx=h.DisplayName; delete(h); end; xyzlabel(sprintf('X (mm) %s',sx),'','','ratio');
        for j=1:8; ha(j).XLim=[-1 1]*500; end;
        for j=[5 6 7 8]; ha(j).YLim=YL(j,:); end;
        for j=[2 4 6 8]; ha(j).YLim(1)=0; end;
        text(ha(2),-1200,ha(2).YLim(2)*1.14,f);  fcnlinewidth(1.5);  fcnfontsize(17);
        export_fig('-a1','-q90','-r250',[A.filename '.png'])
        
        fig(1,1,1.5); for j=1:2; X=inputs(:,vi{j})'; X(257:1024,:)=0;  fcnpulsetemplate(X(:,i),t(i,j)); end; fcntight; xyzlabel('T (sample)','Amplitude (normalized)','',f)
        h=get(gca,'Children'); legend(sprintf('CH1 LEFT%s',h(1).DisplayName),sprintf('CH2 RIGHT%s',h(2).DisplayName),'Location','best'); legend boxoff
        %export_fig('-a1','-q90','-r300',[f '.template.png'])
    end
    
    
    
    
    
end
end


function []=exp1fit(h,xa,b)
F=fit(double(xa),double(b),'exp1');  Fb=1./abs(F.b)/10;
h.DisplayName=sprintf('%.0f mV, %.0f cm',mean(b)*1E3,Fb);
end


function []=exp2fit(h,xa,b)
F=fit(double(xa),double(b),'exp2');  Fb=1./abs(F.b)/10;  Fd=1./abs(F.d)/10;
if Fd>Fb 
    h.DisplayName=sprintf('%.0f mV, %.0f-%.0f cm',mean(b)*1E3,Fb,Fd);
else
    h.DisplayName=sprintf('%.0f mV, %.0f-%.0f cm',mean(b)*1E3,Fd,Fb);
end
end