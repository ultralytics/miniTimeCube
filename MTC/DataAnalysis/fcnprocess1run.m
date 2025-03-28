% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [A, Cout] = fcnprocess1run(A,handles,plotflag)
try closeallexcept(handles.GUI.figure1); end %#ok<TRYNC>
if ~exist('plotflag','var') || isempty(plotflag); plotflag=false; end
A.SRCCHi = MTCpixelID2SRCCH(1:1536);  A.uhj=load('MTCmapping.mat'); A.uhj=A.uhj.X;
pid = int16(MTCSRCCH2pixelID(A.x(:,3:6)));  nw=512;  if eval(A.filename(5:8))==99; nw=16; end
RW=A.x(:,7);  W=A.x(:,8);  i=W<RW;  Wa=W-RW+1;  Wa(i)=Wa(i)+nw;  clear RW

[~, i] = fcnpruninglist(pid);  i=true(size(pid));
%i = i & Wa>416 & Wa<424;  %fig; histogram(Wa,400:512); set(gca,'YScale','log')

% i = i ...
%     & ~(A.SRCCHi(pid,12)==6 & A.x(:,1)==1) ...
%     & ~(A.SRCCHi(pid,12)==3 & A.x(:,1)==2) ...
%     & ~(A.SRCCHi(pid,12)==4 & A.x(:,1)==3) ...
%     & ~(A.SRCCHi(pid,12)==2 & A.x(:,1)==4) ...
%     & ~(A.SRCCHi(pid,12)==1 & A.x(:,1)==5) ...
%     & ~(A.SRCCHi(pid,12)==5 & A.x(:,1)==6); %PRUNE ONLY 1 LASER FACE (UH123456 = UL612435)

pid = pid(i);  if isempty(pid); C=[]; return; end
A.x=A.x(i,:);  W=W(i);  V=A.x(:,10:73)';  Wa=Wa(i);  A.E=A.E(i);  A.unixTime=A.unixTime(i);   A.triggers=A.x(:,9); % Window, nx64 Voltages, Scrod, WindowAdjusted, Event
V = subtractpedestals(A.pedestals,A.pedestalOutliers,pid,W,V);  
A.pedestals=[]; A.pedestalOutliers=[];  clear W

%V=single(V);  fig; fhistogram(V,-80:80); xyzlabel('values (samples)','','',sprintf('%s with 2668 (HVOFF) Pedestal',str_(A.filename)));

%i = ~all(V==1600);  pid=pid(i);  V=V(:,i);  Wa=Wa(i);  A.E=A.E(i);  A.unixTime=A.unixTime(i);  A.triggers=A.x(i,9);  %A.x=[];
%V=removeSpikes(single(V),50);

% d = fdesign.lowpass('Fp,Fst',0.001,.15);
% H = design(d, 'equiripple');  % Design an FIR equiripple filter
% H.Arithmetic = class(single(1));

A.nw=double(max(Wa)-min(Wa)+1); %window adjusted, fwi=first window index
A.rows=numel(A.E);  [A.events, A.ei] = fcnunique(A.E);  A.ei=[A.ei; A.rows+1];  ne=numel(A.events);  clear E %event start indices
V=V';

ne=min(90000,ne);  A.fwi=zeros(ne,1);  C=cell(ne,2,2); C(:,1)={sparse(A.nw*64,1536)}; C(:,2)={sparse(1536,6)};  tic
for i=1:ne
    fprintf('\n%.0f/%.0f  ',i,ne)
    j = A.ei(i):A.ei(i+1)-1;
            
%     if (max(Wa(j))-min(Wa(j)))>7 %AT LEAST 4 WINDOWS
%        j=firstfour(Wa(j),j); %FIRST 4
%        %j=lastfour(Wa(j),j); %LAST 4
%     else
%         j=firstfour(Wa(j),j);
%         %continue %ONLY 1 BOOKEND, AMBIGUOUS
%     end

    Wj=Wa(j);  A.fwi(i)=min(Wj);  nw=max(Wj)-min(Wj)+1;  nj=numel(j);
    
    %[c,~,indC,wpc] = fcnunique(pid(j)); %wpc = windows per channel
    %gc = sparse(double(Wj),indC,1,A.nw,numel(c));  %gap channels (missing window between 2 good windows)
    %s=sum(diff(gc,1,1)==1)'<2 & wpc>1 & wpc<16;  j=j(s(indC));  if isempty(j); continue; end
    
    xa = zeros(nw*64,1536,'int16');
    k = sub2ind([64 nw 1536],ones(nj,1),Wj-min(Wj)+1,pid(j));
    xa(k+(0:63)) = V(j,:); %fig; plot(xa)
    tv=accumarray(pid(j),A.triggers(j),[1536,1]); %trigger vector
    %xa=filter(H,single(xa));
    if i<10000
        [C{i,2}, C{i,1}]=geteventstats(A,xa,tv);  %fig; plot(C{i,1})
    else
        C{i,2}=geteventstats(A,xa,tv);
    end
    
    if i<30 %PLOT SINGLES
       % in=evalin('base','input'); in.MTC.A=A;  in.MTC.C=C;  in.eventNumber=i-1;  eventViewer(in,[]); %export_fig(gcf,'-q90','-r200','-a1',sprintf('f%g.pdf',i));
    end
end; fprintf('\n Done. (%.1fs)\n',toc);

%fcngetQEmap(C);
%fcnsolvetimeoffsets_128_12(A,C)

if plotflag;  fcnplotrun(A,C);  end

%a=[]; for i=1:ne; s=C{i,2,1}(:,4); s(s==0)=nan; s=fcnsigmarejection(s); a(i,:)=nanstd(s); end; fig(1,1,'19cm'); histogram(a,linspace(0,150,200)); xyzlabel('time 1\sigma per event (samples)','','',str_(A.filename))

Cout = C;
end

function [stats, x] = geteventstats(A,x,tv)
[j, a, w, t, integral, x, minx] = fcnpulsewidth(single(x),[.45 .50 .55],[100 90000],[14 21],1600,'fraction');  k=t~=0; %[10 128] widths before
%t(k) = mod(t(k)+(A.fwi-1*0)*64,128); %128 vector time

%APPLY CALIBRATIONS
a(j) = a(j)./A.offsets(j,2);                            a(isnan(a) | a<0)=0;
w(j) = w(j)./A.offsets(j,3);                            w(isnan(w) | w<0)=0;
t(k) = t(k) - A.offsets(k,4);                           t(isnan(t))=0;
integral(j) = integral(j)./A.offsets(j,5);              integral(isnan(integral))=0;

if nanstd(fcnsigmarejection(t(t~=0)))>15; t=t*0; w=w*0; end

%fig;  a=dt; a(a==0)=nan; fcnplotdetectorprojection(input,ones(1536,1),a,[0 0 0],'winkeltripel','SCROD'); %colormap('bluewhite')
stats = single([j, a, w, t, integral, tv, minx]); %[triggered, amplitude, width, time, integral, triggerKurtis]
end

function fcnsolvetimeoffsets_128_12(A,C)
z = cell2mat(C(:,2));  z=z(:,4);  i=find(z);
j = mod(i,1536);  j(j==0)=1536;  %j = cell2mat(C(:,3));  j=j(i);  %j=pid 1-1536
z=z(i);  nz=numel(z);
si=A.SRCCHi(j,10);  ci=A.SRCCHi(j,11);  %channel indices on scrod (1-128)

r = repmat((1:nz)',[1 2]); %rows
c = [si, ci+12];
H = sparse(r(:),double(c(:)),1,nz,12+128);

i=sum(H)>0;  H=H(:,i);  x=zeros(12+128,1);
x(i) = (H'*H)\H'*z; %LLS
%x = zeros(12+16+8,1); for i=1:6;  x = x + (H'*H)\H'*(z-x(si)-x(ai+12)-x(ci+28)); end %NLS
% R = speye(nz,nz)*.0384^2; %measurement noise matrix
% C = inv(H'/R*H); %covariance matrix
% s = sqrt(diag(C)); %sigmas

ui=unique(si);  so=zeros(12,1);   so(ui)=x(ui   )-mean(x(ui   ));
ui=unique(ci);  co=zeros(128,1);  co(ui)=x(ui+12)-mean(x(ui+12)); %save MTCDT128map co;
%fig; fcnplotdetectorprojection(evalin('base','input'),1:1536,so(A.SRCCHi(:,10)),[0 0 0],'winkeltripel','SCROD'); h=colorbar('eastoutside','AxisLocation','out'); h.Label.String='Time Offset (bins)'; colormap(redblue); title(sprintf('128 Time Offsets from ''%s''',str_(A.pf2)))
%fig; fcnplotdetectorprojection(evalin('base','input'),1:1536,co(A.SRCCHi(:,11)),[0 0 0],'winkeltripel','SCROD'); h=colorbar('eastoutside','AxisLocation','out'); h.Label.String='Time Offset (bins)'; colormap(redblue); title(sprintf('128 Time Offsets from ''%s''',str_(A.pf2)))
%input = evalin('base','input'); fig; fcnplotdetectorprojection(input,1:1536,input.MTC.A.po,[0 0 0],'winkeltripel','SCROD'); h=colorbar('eastoutside','AxisLocation','out'); h.Label.String='Time Offset (bins)'; colormap(redblue); title(sprintf('1536 7/13/2015 CAJIPCI-Pulse Time Offsets'))
end


function fcngetQEmap(C)
ne = size(C,1);  a = reshape(full(cell2mat(C(:,2))),[1536 ne 7]);
d=mean(a(:,:,1),2); %trigger rate
a(:,:,4) = mod(a(:,:,4),64); %mitigate window slipping effect
a(a==0)=nan;

map = nan(1536,5);  map(:,1)=d;
for i=1:1536
    for j=2:5
        map(i,j) = mean(fcnsigmarejection(a(i,:,j)));
    end
end
map(:,4)=map(:,4)-nanmean(map(:,4));

save -v6 MTCoffsetsNew.mat map
s={'trigger','amplitude','width','time','integral'};
ha=fig(5,1); for i=1:5; h=sca; fcnplotdetectorprojection(input,1:536,map(:,i),[0 0 0],'winkeltripel','PMT'); hc=colorbar('East','AxisLocation','out','Location','southoutside'); hc.Label.String=s{i}; h.CameraViewAngle=4.9053; end; colormap(ha(1),bluewhite)
fig(2,2); for ch=1; b=squeeze(a(:,:,2:5)); for i=1:4; sca; fhistogram(b(:,:,i),200); xyzlabel('','','',s{i+1}); end; end
end


function []=fcnplotrun(A,C)
type = 'PMT';
nc=5; %number of columns
mof = false; %maps only flag

uhj=sortrows(A.uhj,10);  input=evalin('base','input'); 
switch type %13=255, 102=0
    case 'SCROD';           va=unique(uhj(:,2),'stable');   tid = A.SRCCHi(:,10);
    case 'ASIC row';        va = 0:3;    tid = A.SRCCHi(:,2)+1;
    case 'ASIC column';     va = 0:3;    tid = A.SRCCHi(:,3)+1;
    case 'ASIC channel';    va = 0:7;    tid = A.SRCCHi(:,4)+1;
    case {'PMT','MCP'};     va = 1:24;   tid = A.SRCCHi(:,5);
    case 'ASIC';            va = 1:16;   tid = A.SRCCHi(:,8);
    case 'Detector ASIC';   va = 1:192;  tid = A.SRCCHi(:,9);
end;  tid=double(tid);
ne = size(C,1);
na = numel(va);  ov=ones(1536,1);       titlestr=sprintf('%.0f Events from %s',ne,str_(A.filename));

str = {'Triggers','Amplitude','Width','Time','Integral','Candidate','Minima'};
stru = {'','samples','samples','samples','','',''};
sf = {' %0.2f',' %.0f',' %.1f',' %.1f',' %.0f',' %.2f',' %.1f'};
a = reshape(full(cell2mat(C(:,2))),[1536 ne 7]);  pl = fcnpruninglist;
izf = [1 0 0 0 0 1 0]==1; %ignore zeros flag
if mof; ha=fig(1,5,'10x60cm');  else ha=fig(3,nc,'22.5x80cm');  end %#ok<UNRCH>
for ip = 1:nc
    b = double(a(:,:,ip));
    x = fcnaccumrows(b,tid,ov,[na ne]);  if ~izf(ip); d=b~=0; else d=ones(size(b)); end
    n = fcnaccumrows(d,tid,ov,[na ne]);  n(n==0)=1;  x=x./n;
    x(x==0)=nan;  c=nanmean(x,2);  c(c==0)=nan;

    if mof
        sca(ha(ip));  e=sum(b,2)./sum(d,2); if ~izf(ip); e(e==0)=nan; end
        fcnplotdetectorprojection(input,pl,e,[0 0 0],'winkeltripel',type); hc=colorbar('East','AxisLocation','out','Location','southoutside'); hc.Label.String=sprintf('%s',stru{ip}); xyzlabel(gca,'','','',sprintf('%s',str{ip}));
        if izf(ip); colormap(gca,bluewhite); end
        if ip==3; title(sprintf('%s - %s',str_(A.filename),str{ip})); end
        continue
    end
    %y=b(b~=0); fig; histogram(y,200); fcntight; set(gca,'YScale','log'); %DARK COUNTS
    
    sca(ha(ip));  [su, sl]=fcnstd(x');  fcnerrorbar(1:na,c,su,sl,'b');
    if na<25; for i=1:na; text(i,c(i)+su(i),sprintf(sf{ip},c(i)),'horizontalalignment','left','verticalalignment','middle','rotation',90); end; end
    xyzlabel(type,stru{ip},'',sprintf('%s',str{ip})); if ip==3; title(sprintf('%s\n%s',str_(A.filename),str{ip})); end
    sca(ha(ip+nc)); pcolor((1:na+1)-.5,1:ne,x([1:na na],:)'); shading flat; xyzlabel(type,'event'); grid off; %colorbar('Location','East','color',[1 1 1]*0);
    fprintf('Average %s: %.2f %s\n',str{ip},nanmean(c),stru{ip})
       
    %am=a(:,:,2);
    %w=a(:,:,3);
    %an=a(:,:,7);
    %i=am>500 & w~=0;
    %hist211(w(i),am(i),200); xyzlabel('Width (bins)','Amplitude (PE)')
    %hist211(am(i),an(i),{linspace(0,1400,200),linspace(-200,0,200)},'scatter'); xyzlabel('Max (bins)','Min (bins)')
    %hist211(am(i),w(i),{linspace(0,1400,200),linspace(0,30,200)}); xyzlabel('Amplitude (bins)','Width (bins)')

    e = zeros(1536,1);
    if ~izf(ip)
        for j=1:1536
            es=b(j,:); e(j)=mean(fcnsigmarejection(es(es~=0),3,3));
            %e(j)=median(es(es~=0));
        end
    else
        e=sum(b,2)./sum(d,2);
    end
    sca(ha(ip+nc*2));                           if ~izf(ip); e(e==0)=nan; end    


    %sca(ha(ip+nc*2));  e=sum(b,2)./sum(d,2);   if ~izf(ip); e(e==0)=nan; end
    fcnplotdetectorprojection(input,pl,e,[0 0 0],'winkeltripel',type); colorbar('eastoutside','AxisLocation','out'); set(gca,'CameraViewAngle',4.6);
    ha(ip+nc).CLim = ha(ip+nc*2).CLim;
    map(:,ip)=e;

    
    ribbonFlag=false;
    if ribbonFlag && ip==2
        xi=linspace(0,200,50);
        yh=zeros(na,50);
        for i=1:na
            y=b(tid==va(i),:); y=y(y~=0);
            [yhi, xh] = fcnhistc(y,xi);
            yh(i,:) = yhi;
        end
        fig; fcnribbon(xh(1:end-1),yh(:,1:end-1)); xyzlabel('PMT','V (bins)','',titlestr);
    end

    if izf(ip); colormap(ha(ip+nc),bluewhite); colormap(ha(ip+nc*2),bluewhite); end
end
%map(:,4)=map(:,4)-nanmean(map(:,4));  mapzeros=map; mapzeros(isnan(mapzeros))=0;  save MTCoffsets.mat map mapzeros
%fig(1,5); for i=1:5; sca; [y,x]=fcnhistc(a(:,i),50,fcndefaultcolors(1)); histNfit(x,y,'b'); xyzlabel(stru{i},'','',str{i}); end %6.27 and 6.30 difference histograms

if ~mof;  fcntight(ha(1:nc*2),'jointx');  fcntight(ha(nc:nc*2),'y');  if na<25; set(ha,'xtick',1:na,'xticklabel',va); end; end
if mof
    text(.5,1.2,titlestr,'FontSize',14,'HorizontalAlignment','left','Units','Normalized','Parent',ha(1)); 
    h=findobj(gcf,'type','colorbar'); for i=1:numel(h); h(i).Position(4)=.03;  h(i).Position(2)=.15; h(i).Box='off'; end
    set(ha,'CameraViewAngle',5.8); fcnfontsize(16)
end

end

function combineWaveforms()
a=zeros(1,1536); X=zeros(size(input.MTC.C{1,1}));
for i=1:size(input.MTC.C,1)
    b=input.MTC.C{i,1};  c=input.MTC.C{i,2};  
    %j=find(c(:,4)~=0);  b=b(:,j); c=c(j,:);
    b=b./max(max(b),1);
    X=X+b;  a=a+any(b,1);
end
end



function SPR()
n=size(input.MTC.C,1);  X=cell(n,1);
for i=1:n
    b=double(input.MTC.C{i,1});  c=input.MTC.C{i,2};  
%    j=find(c(:,4)==0);  b=b(:,j); c=c(j,:);  t=1:size(b,1);
%     for j=1:size(c,1)
%         b(:,j)=interp1c(t,b(:,j),t+c(j,4)-150);
%     end

    if max3(b)<100; b(:)=nan; end
    X{i}=b;
    i
end
X=cell2mat(X); X=reshape(X,[256 n 1536]); X=permute(X,[3 1 2]);  mu=nanmean(X,3);  s=nanstd(X,[],3);
%ha=fig; pcolor(mu); shading flat; fcntight; ha.CLim=[-30 30];  save laserSPRTemplate.mat mu

[~, ~, ~, ~, PMT, PMTR, PMTC] = MTCpixelID2SRCCH(1:1536);


close all
border = [0, 0, 1.9, 0.6, 1.6, 1.2];
%border = [spacingHoriztontal, spacingVertical, leftBorder, rightBorder, bottomBorder, topBorder];
ha = fig(8,8,.4,.4,border);  [~,brightest]=max(max(mu,[],2));
for r=1:8
    for c=1:8
        h=sca(ha((r-1)*8+c));
        j=PMTR==r & PMTC==c & PMT==21;
        plot(squeeze(X(j,:,:)))
        if j(brightest); color=fcndefaultcolors(2); else; color=fcndefaultcolors(1); end
        %he=errorarea(1:256,mu(j,:),s(j,:),color); 
        grid off; box on;
        if c>1; h.YTick=[]; h.YColor='k'; else; h.YTick=[0 1600 round(max3(mu))]; end
        if r<8; h.XTick=[]; h.XColor='k'; else; h.XTick=[256]; end
    end
end
fcntight('x'); fcntight('yjoint'); fcnlinewidth(1)
ha(4).Title.String = sprintf('%g events from ''%s''',n,str_(input.MTC.A.filename));
%for i=1:64; ha(i).YLim(2)=1600; end
end


%FIND OUTLIERS IN Ts and Qs
function []=fcnTQ()
s=[]; mu=[]; ne=numel(input.MTC.events); 
for i=1:ne 
    a=input.MTC.C{i,2};  k=find(a(:,4)~=0);
    [~,j]=fcnsigmarejection(a(k,4),3,2);  a(k(~j),4)=0; %reject bad times

    a(:,1)=sum(a(:,4)~=0);
    a(:,6)=sum(a(:,6)~=0);
    b=a(a(:,4)~=0,:);
    s(i,:)=nanstd(b,[],1);  
    mu(i,:)=nanmean(b,1); 
end
str={'Triggers','Amplitude','Width','Time','Integral'}; ha=fig(2,5);
j=true(size(mu,1),1); for i=1:5; [~,jmu]=fcnsigmarejection(mu(:,i),6,3);  [~,js]=fcnsigmarejection(s(:,i),5,3); j=j & jmu & js;  hv{i,1}=linspace(0,max(mu(j,i)),50);  hv{i,2}=linspace(0,max(s(j,i)),50);end
for i=1:5; sca(ha(i)); histogram(mu(:,i),hv{i,1}); xlabel([str{i} ' \mu']); axis tight; sca(ha(i+5)); histogram(s(:,i),hv{i,2}); xlabel([str{i} ' \sigma']); axis tight; end
title(ha(3),sprintf('%s %g DELAYED events',str_(input.MTC.A.filename),ne))
%save PROMPT.mat mu s

clear a b
load PROMPT.mat; a.mu=mu; a.s=s;
load DELAYED.mat; b.mu=mu; b.s=s;
hist211(b.s(:,4),a.mu(:,2),50); xyzlabel('DELAYED Time \sigma','PROMPT Energy (PE/pixel)')
end