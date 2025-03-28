% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function []=eventViewer(input,output)
type                    = 'PMT';
meanflag                = false;
normalizeflag           = false;
applyCalibrationsFlag   = false;

handles = evalin('base','handles');
if nargin==0
    input = evalin('base','input');
    output = evalin('base','output');
end
MTCflag = ischecked(handles.GUI.realdataflag);
LAPPDflag = evalin('base','flags.status.lappd');

hp=findobj('type','figure','Name','Event Viewer');
if ~isempty(hp)
    hf=hp;  ha=flip(hf.Children(isaxis(hf.Children)));  for i=1:numel(ha); cla(ha(i)); end
else
    [ha,hf]=fig(2,3,1.2,1.6);  hf.Name='Event Viewer';
    ha(1).Position(2)=ha(4).Position(2); ha(1).Position(4)=.88;  delete(ha(4));  ha=ha([1 2 3 5 6]);
    if MTCflag; xyzlabel(ha(1),'T (Window)','PMT'); xyzlabel(ha(3),'T (Window)','V (bins)');
    else xyzlabel(ha(1),'T (bin)','Channel'); xyzlabel(ha(3),'T (ns)','V');  end
    xyzlabel(ha(2),'','','','Amplitude',16);
    xyzlabel(ha(4),'','','','Time',16);
    xyzlabel(ha(5),'','','','Integral',16);
    h=colorbar(ha(1),'Axislocation','out');  h.Title.String='';  h.Title.VerticalAlignment = 'bottom'; 
    h=colorbar(ha(2),'Axislocation','out');  h.Title.String='';  h.Title.VerticalAlignment = 'bottom';
    h=colorbar(ha(4),'Axislocation','out');  h.Title.String='';  h.Title.VerticalAlignment = 'bottom';  %colormap(ha(4),'redblue')
    h=colorbar(ha(5),'Axislocation','out');  h.Title.String='';  h.Title.VerticalAlignment = 'bottom';
end

h=ha(1);
if MTCflag
    A=input.MTC.A;  C=input.MTC.C;  ei=input.eventNumber+1;
    y=C{ei,1};   %y=y((257-257)+(1:256),:);
    s=full(C{ei,2});  SI=sortrows(A.uhj,[2 6]);  
    %[smu,ssig]=weightedMean(s,repmat(s(:,4)~=0,[1 6])) %[triggers amplitude width time integral]
    %i=s(:,4)~=0; s(i,:); fig; fcnhist(s(i,4))
    %std(fcnsigmarejection(s(i,4),3,3))
    
    %v=A.SRCCHi(:,5)~=21;  y(:,v)=0;  s(v,:)=0; %PMT 21 ONLY
    %v=any(A.SRCCHi(:,1:3)~=[4 1 2],2);  y(:,v)=0;  s(v,:)=0; %SRC 4-1-2 only
    %v=~(s(:,4)>0 & s(:,6)==0);  y(:,v)=0;  s(v,:)=0; %WAVEFORMS NO TRIGGERS
    %v=~(s(:,4)==0 & s(:,6)>0);  y(:,v)=0;  s(v,:)=0; %TRIGGERS NO WAVEFORMS
    
    x=find(any(y,2));
    if any(x); x=x(1):x(end); end;  
    y=full(y(x,:));   
    y(y==0 & circshift(y,1,1)==0 & circshift(y,-1,1)==0)=nan;  %laserTemplate=load('laserTemplate.mat');  %y=y-laserTemplate.mu';
    y0=y;  y=y(:,SI(:,1)); x=x+A.fwi(ei,1)*64;  if isempty(y); minx=0; maxx=1; else minx=min(x); maxx=max(x); end

    if any(x); pcolor(h,x,1:1536,y'); shading(h,'flat'); end;
    h.CLim=[-30 30]; %[-500 1000];
    xyzlabel(ha(3),'','','',sprintf('%s event %.0f PROMPT',str_(A.filename),A.events(ei)));  
    mmx = fcntight(h,'xy');
    
    %vi = unique(SI(:,2),'stable');  yi=mean(reshape(1:1536,[128 12]))-63;   h.YTick=yi; h.YTickLabel=vi;
    vi = unique(SI(:,6),'stable');  yi=mean(reshape(1:1536,[64 24]))-31.5;   h.YTick=yi; h.YTickLabel=vi;
    vi = round(mmx/64); vi=vi(1):vi(end);  h.XTick=vi*64;  h.XTickLabel=vi;  
    n=numel(vi); if n>12; i=true(n,1); i(round(linspace(1,n,12)))=false; h.XTickLabel(i,:)=' '; end
    
    
    if applyCalibrationsFlag
        s(:,2)=s(:,2)./A.offsets(:,2); %Amplitude Map
        i=s(:,4)~=0;  s(i,4)=s(i,4)-A.offsets(i,4);  %Time Map
        s(:,5)=s(:,5)./A.offsets(:,5); %Integral Map
    end
    
    h=ha(3);
    SRCCHi = A.SRCCHi; %MTCpixelID2SRCCH(1:1536);
    switch type %13=255, 102=0
        case 'SCROD';           uhj=sortrows(A.uhj,10);  va=unique(uhj(:,2),'stable');   tid = SRCCHi(:,1);
        case 'ASIC row';        va = 0:3;    tid = SRCCHi(:,2);
        case 'ASIC column';     va = 0:3;    tid = SRCCHi(:,3);
        case 'ASIC channel';    va = 0:7;    tid = SRCCHi(:,4);
        case {'PMT','MCP'};     va = 1:24;   tid = SRCCHi(:,5);
        case 'ASIC';            va = 1:16;   tid = SRCCHi(:,8);
        case 'Detector ASIC';   va = 1:192;  tid = SRCCHi(:,9);
    end;
    na = numel(va);  c=fcndefaultcolors(1:na,na);  ap=any(y0,1)'; %active pixels
    for i=1:na
        ci = tid==va(i) & ap; %channel indices
        yi = y0(:,ci);  nc=size(yi,2);
        xi=repmat(x',[1 nc]);
        if any(isfinite(yi(:)))
            if applyCalibrationsFlag
                ti = xi - A.offsets(ci,4)'; 
                if meanflag
                    for j=1:nc
                        yi(:,j) = interp1(ti(:,j),yi(:,j),xi(:,j),'linear',0);
                    end
                else
                    xi=ti;
                end
            end
            
            if meanflag;        yi=nanmean(yi'); xi=x;              else  xi(end+1,:)=nan; yi(end+1,:)=nan; end
            if normalizeflag;   yi=yi./max(yi);     h.YLabel.String='Normalized';   end
            plot3(h,xi(:),yi(:),zeros(numel(yi),1)+i,'-','Color',c(i,:),'DisplayName', sprintf('%s %g',type,va(i)),'linewidth',1.5);
        end
    end
    h.XTick=ha(1).XTick; h.XTickLabel=ha(1).XTickLabel;  h.XLim = [minx maxx];  %if ~normalizeflag; h.YLim=[-1000 1500]; else h.YLim=[-.5 1];  end
    %plot(wbx,wby,'-','color',[1 1 1]*.3)
else
    D=output(1).D;  
    if LAPPDflag %LAPPD
        sd=size(D); input.cube.pixels = prod(sd(2:end));
        D=reshape(D,[sd(1) input.cube.pixels]);
    else
        s = output(1).ZI;
        %E=D; E(E==0)=nan;  %D=floorandceil(D,-1.1,1.1);
        %s(:,5)=nanvar(E,[],1); %VARIANCE
       % s(:,2) = sum(abs(D),1); %ABSOLUTE INTEGRAL
        %s(:,2) = rms(D,1); %RMS
    end
    pcolor(h,1:input.cube.dsp.samples,1:input.cube.pixels,D'); shading(h,'flat');

    D(end+1,:)=nan;  x=input.cube.dsp.t;  x(end+1)=nan;  x=repmat(x',[1 input.cube.pixels]);
    plot(ha(3),x(:),D(:),'.-','Markersize',.1); 
    fcntight(ha([1 3]),'xy')
    type='';
    xyzlabel(ha(3),'','','',sprintf('%s event %.0f',str_(evalin('base','GEANT.pf1')),input.eventNumber));  
end

if ~LAPPDflag
    sca(ha(2));  a=s(:,2);  a(a==0)=nan;  fcnplotdetectorprojection(input,find(a),a,[0 0 0],'winkeltripel','');  b=fcnsigmarejection(a,3,3);  ha(2).Title.String=sprintf('Amplitude \n%.2g \\pm %.2g',mean(b),std(b));
    sca(ha(4));  a=s(:,4);  a(a==0)=nan;  fcnplotdetectorprojection(input,find(a),a,[0 0 0],'winkeltripel','');  b=fcnsigmarejection(a,3,3);  ha(4).Title.String=sprintf('Time\n%.3g \\pm %.2g',mean(b),std(b));   %fig; fcnhist(rem(b-median(b),64))
    sca(ha(5));  a=s(:,3);  a(a==0)=nan;  fcnplotdetectorprojection(input,find(a),a,[0 0 0],'winkeltripel','');  b=fcnsigmarejection(a,3,3);  ha(5).Title.String=sprintf('Width\n%.2g \\pm %.2g',mean(b),std(b));
    set(ha([2 4 5]),'CameraViewAngle',6.5);  
    linkaxes(ha([2 4 5]));
    %ha(1).CLim(2) = ha(2).CLim(2);
    
    %fig; A=input.MTC.A; a=double(A.SRCCHi(:,12)==6); fcnplotdetectorprojection(input,ones(1536,1),a,[0 0 0],'winkeltripel','SCROD'); colormap('bluewhite')
    %fig;  load('MTCoffsets.mat'); a=map(:,1); a(a==0)=nan; fcnplotdetectorprojection(input,ones(1536,1),a,[0 0 0],'winkeltripel','PMT'); %colormap('bluewhite')
    
%     fig; [~, a] = fcnpruninglist(); a=double(a); a(a==0)=nan; fcnplotdetectorprojection(input,ones(1536,1),a,[0 0 0],'winkeltripel','PMT'); %colormap('bluewhite')
    
    
    %NTC SCROD
    %fig;  a=ones(128,1); a(a==0)=nan; fcnplotdetectorprojection(input,ones(128,1),a,[0 0 0],'winkeltripel','channel'); %colormap('bluewhite')
end