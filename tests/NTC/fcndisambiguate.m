function results = fcndisambiguate(input,output,handles,PE,flags,plotflag)
if plotflag; closeallexcept(handles.GUI.figure1);  deleteh(findobj(handles.GUI.figure1,'DisplayName','test')); end
A=output(1).D;
[nx, ns, ~, nlappd] = size(A); %number of times, number of strips (30), number sides (2), number of lappd's

Xa = repmat((1:nx)',[1 ns]);
Ya = repmat(1:ns,[nx 1]);
bi=1:1:256;  nx=numel(bi);  Xa=Xa(bi,:); Ya=Ya(bi,:); A=A(bi,:,:,:);  input.cube.dsp.dt=0.1;  A0=A;

xb=1:.1:nx;
yb=1:.1:ns;
[Xb,Yb] = ndgrid(xb,yb);

%PULSE RECONSTRUCTION
B=reshape(A,[nx ns*2*nlappd]);
B=fcnSplinePulseReconstruction(B,input.cube.dsp.dynamicrange(2));  Ad=reshape(B,[nx ns 2 nlappd]);

%DECONVOLVE DOWN STRIPS
nsr = .004^2/var(B(:));
[psf,S]=fcngetalongpsf(input,flags);  B=deconvwnr(B',psf,nsr);  A=reshape(B',[nx ns 2 nlappd]);
S=deconvwnr(S,psf,nsr); dSPRamplitude=max(S); %deconvolved SPR amplitude
%threshold = max([ max(B(:))*0.2 .05]); %Volts
threshold = max([ max(B(:))*0.05 .005]); %Volts

% hb=fig(3,1,'28.5x19cm'); %3 plots for paper
% sca; plot3(Xa,Ya,A0(:,:,1,1),'.-'); %Signal
% sca; plot3(Xa((1:64)+96,(1:11)+7),Ya((1:64)+96,(1:11)+7),psf2,'.-');  %PSF
% %sca; surf(hb(3),Xb,Yb,Zb); shading flat;  alpha(.8);
% fcnlinewidth(1.5); fcnmarkersize(15); fcnfontsize(16); fcnview(hb,'skew'); fcntight(hb,'xyz'); set(hb(2),'Xlim',[1 256],'Ylim',[1 26]); set(hb,'XDir','reverse','Ydir','reverse'); xyzlabel(hb,'T (sample)','Strip (1-26)','Amplitude (V)')

%DECONVOLVE ACROSS STRIPS
B = reshape(permute(A,[2 3 4 1]),[ns 2*nx*nlappd]);
psf=fcngetacrosspsf(input);  B=deconvwnr(B',psf,nsr);  A=permute(reshape(B',[ns 2 nlappd nx]),[4 1 2 3]);

%fig; x = linspace(-100,100,26); xb=linspace(-100,100,200); ya=fcnpdfnorm(x,-47,5);  plot(xb,fcnpdfnorm(xb,-47,5),'-','color',[1 1 1]*.7);
%F = griddedInterpolant(x,ya,'cubic');  yb=F(xb); plot(x,ya,'k.','Markersize',30); plot(xb,yb,'r.')

for tile=1
    if plotflag
        ha=fig(2,2,'19x19cm');        
        sca(3); plot3(Xa,Ya,A0(:,:,1,tile),'.-');  title('LEFT')
        sca(4); plot3(Xa,Ya,A0(:,:,2,tile),'.-');  fcnview(ha(3:4),'skew'); set(ha(3:4),'XDir','reverse','Ydir','reverse'); title('RIGHT')
    end
    for side=1:2
        Za = A(:,:,side,tile);
        F=griddedInterpolant(Xa,Ya,Za,'spline');  Zb=F(Xb,Yb);
        [px,py,z]=fcnfind2Dpeaks(Zb,threshold); 

        D=[xb(px)' yb(py)' z];  if plotflag; sca; pcolor(Xb,Yb,Zb); shading flat;  alpha(.7);  plot(xb(px),yb(py),'k+'); end
        if side==1; L=D; else R=D; end
    end
    Ps=fcnPs(L,R); %strip probability
    Pa=fcnPa(L,R); %amplitude probability
    Pt=fcnPt(L,R); %time probability
    P = Ps.*Pa.*Pt;
    [nl,nr]=size(P);  if nl==0 || nr==0 || max3(A(:,:,1,tile))<.025; fprintf('No peaks found on LAPPD%g\n',tile); continue; end
    
    %MATCH
    [L,R]=matchLR(L,R,nl,nr,P);
    n = size(L,1);  if n==0;  fprintf('No matches found on LAPPD%g\n',tile); continue; end
    
    if plotflag
        c=[1 1 1]*1;
        sca(1); ylabel('Strip (1-26)'); plot(L(:,1),L(:,2),'ko');  for i=1:n;  text(L(i,1),L(i,2),sprintf(' %.0f',i));  fcnfontsize(gca,14);  end; text(10,22,'LEFT','fontsize',35,'color',c); 
        sca(2); ylabel('Strip (1-26)'); plot(R(:,1),R(:,2),'ko');  for i=1:n;  text(R(i,1),R(i,2),sprintf(' %.0f',i));  fcnfontsize(gca,14);  end; text(10,22,'RIGHT','fontsize',35,'color',c)
        fcnmarkersize(ha,8); fcnlinewidth(ha,1);  xyzlabel(ha,'T (sample)','Strip (1-26)','Amplitude (V)');  fcntight(ha,'xyz');  fcntight(ha(1:2),'joint c0');  
    end
    [x,y,t]=fcnLR2xyt(input,L,R,dSPRamplitude);  t = t - (input.cube.pmt.spr*3) - .16; %account for transit time
    xyz=lappd2detector(input,tile,x,y,plotflag);  n=numel(t);
    
    %GET ERRORS
    i=PE.t<25.6 & PE.pixel==tile;  E=[PE.x(i,:) PE.t(i)];  nPE=sum(i);
    r=sqrt((xyz(:,1)'-E(:,1)).^2 + (xyz(:,2)'-E(:,2)).^2 + (xyz(:,3)'-E(:,3)).^2 + (180*(t'-E(:,4))).^2);
    [~,i] = min(r,[],1);
    error = [xyz t]-E(i,:); if n>1; sigmas=std(error,0,1); else sigmas=error; end
    %fig; [y,x]=fcnhist(error(:,4),50,'b','bar');  histNfit(x,y,'b'); xlabel('error (ns)')
    
    %fig; fcnhist(t-E(i,4),30)
    results.xhat = [n sigmas];
    results.true = [nPE 0 0 0 0];
    results.PEhat = [xyz t];
    results.PEtrue = E;
    pr=mean(sigmas(2:3)); tr=sigmas(4); 
    fprintf('LAPPD%g: %g matches (%gL - %gR) from %gPEs (%.1fmm, %.3fns, FI=[%.2f %.0f])\n',tile,n,nl,nr,nPE,pr,tr,n/pr.^2,n/tr.^2)
end
end


function []=TSplots()
n = numel(tsv);  ha=fig(n,3,.6,1);
for i=1:n
    a=MC.xhat(:,:,i); 
    j1=a(:,1)==1;   n1=sum(j1);   b1=std(a(j1,:));   if n1==0; b1=zeros(1,5); end
    j2=a(:,1)>1;    n2=sum(j2);   b2=mean(a(j2,:));  if n2==0; b2=zeros(1,5); end;    b=(b1*n1+b2*n2)/(n1+n2);
    fprintf('\\hline \\bf %6.0f & %6.2f	& %6.1f & %6.3f \\\\ \n',[tsv(i) mean(a(:,1))/tsv(i) mean(b(3:4)) norm([b(5) .060])]); 
    
    
    a = cell2mat(MC.PEhat(:,i));
    b = cell2mat(MC.PEtrue(:,i)); m=201;
    xp = midspace(-100,100,m);
    c=fcnsigmarejection(a(:,4),6,1); xt = linspace(min3(c)-.1,max3(c)+.1,50);  eff=numel(a)/numel(b);
    a(:,4) = a(:,4) + randcdfc(input.cube.pmt.transit.c-.180,numel(a(:,4))); %add TTS
    
    sca; [y,~,h]=fcnhist(a(:,2),xp,fcnlightencolor('r'),'line','est'); h.YData=y/eff;  fcnhist(b(:,2),xp,'r','line','true'); xlabel('Along LAPPD (mm)'); ylabel(sprintf('%g PEs',tsv(i)),'fontsize',30)
    sca; [y,~,h]=fcnhist(a(:,3),xp,fcnlightencolor('g'),'line','est'); h.YData=y/eff;  fcnhist(b(:,3),xp,'g','line','true'); xlabel('Across LAPPD (mm)')
    sca; [y,~,h]=fcnhist(a(:,4),xt,fcnlightencolor('b'),'line','est'); h.YData=y/eff;  fcnhist(b(:,4),xt,'b','line','true'); xlabel('t (ns)'); fcntight(gca,'x')
end; fcnlinewidth(1); fcnmarkersize(10); fcntight; title(ha(2),'Scintillation','fontsize',25);
end


function [L,R]=matchLR(L,R,nl,nr,P)
if nl>=nr; P=P'; end
[v,j]=max(P,[],1);
[v,k]=sort(v,2,'descend');  i=v>1E-16; %eps
li = k(i); %left indices
ri = j(k(i)); %right indices

a=[li(:) ri(:)];  [~,l]=fcnunique(ri(:)); a=a(l,:);
if nl<nr; a(:,[1 2])=a(:,[2 1]); end;  L=L(a(:,1),:);  R=R(a(:,2),:);
end


function [x,y,t]=fcnLR2xyt(input,L,R,dSPRamplitude)
pmt = input.cube.pmt;   l=input.cube.pixelsize(1); %length of LAPPD, 200mm or 800mm;
speed = 190; %(mm/ns) pulse speed down anode
%a=solve('(x-0)/speed = (t1-t)','(l-x)/speed = (t2-t)','x','t');
bins2ns = 1/input.cube.dsp.rate;  t1=L(:,1)*bins2ns;  t2=R(:,1)*bins2ns; %input.cube.dsp.rate=10 Gs/s
x = speed*(t1-t2)/2; %mm from center
t = (t1+t2 - l/speed)/2; %ns
y = interp1c(1:pmt.strips,pmt.stripx, (L(:,2)+R(:,2))/2); %mm from center

%BREAK DOWN LARGE PULSES INTO SMALLER SINLGES
a=floor( (L(:,3)+R(:,3))/2/dSPRamplitude ); ma=max(a);
if ma>1
    c=cell(ma,1); tb=c; xb=c; yb=c;  tb{1}=t; xb{1}=x; yb{1}=y;
    for i=2:ma
        j=a==i;
        b=repmat(t(j),[1 i-1]);  tb{i}=b(:);
        b=repmat(x(j),[1 i-1]);  xb{i}=b(:);
        b=repmat(y(j),[1 i-1]);  yb{i}=b(:);
    end
    x=cat(1,xb{:});  y=cat(1,yb{:});  t=cat(1,tb{:});
end
end


function X=lappd2detector(input,pid,x,y,plotflag)
z=y; y=x; x=zeros(size(x)); nx=numel(x);
rpy = repmat(fcnVEC2RPY(-input.cube.all.normalVec(pid,:)),[nx 1]);
X = rotateW2Bc(rpy,[x(:) y(:) z(:)]);
X = X + input.cube.all.xyz(pid,:);

% if plotflag
%     h=evalin('base','handles.GUI.axes1');
%     photons = evalin('base','photons');
%     i=PE.t<25 & PE.pixel==pid;
%     plot3(h,photons.endPos(i,1),photons.endPos(i,2),photons.endPos(i,3),'g.','Markersize',15,'DisplayName','test')
%     plot3(h,X(:,1),X(:,2),X(:,3),'r.','Markersize',15,'DisplayName','test')
% end
end

function P=fcnPs(L,R) %strip probability
sigma = .1; %(delta-strip units)
P = fcnpdfnorm(L(:,2)-R(:,2)',0,sigma);
end

function P=fcnPa(L,R) %amplitude probability
sigma = std([L(:,3); R(:,3)])*.2; %(V)
P = fcnpdfnorm(L(:,3)-R(:,3)',0,sigma);
end

function P=fcnPt(L,R) %time probability
s = .005; %ns
Gss = 10; %Gs/s
mu = 10.526/Gss; %time bins it takes charge to cross 200mm @190mm/ns
dt = abs(L(:,1)-R(:,1)')/Gss;
P = 1-cdf('norm',dt,mu,s);
%P = 1-fcncdfnorm(dt,mu,s);
%x=linspace(0,15,1000)/Gss;  P=1-cdf('norm',x,mu,s);  fig(1,1,'9.5x19cm'); plot(x,P); xyzlabel('dt = abs(t_L - t_R) (ns)','CDF'); fcnlinewidth(3); legend('P_t'); fcnfontsize(16)
end

function [psf,S]=fcngetalongpsf(input,flags)
p=[];
p.t = input.cube.dsp.dt*128;
p.pixel = 1;
p.xpixel = [0 0 input.cube.pmt.stripx(1)];
S = fcnanalogvoltage(input,flags,p,input.cube.dsp.t,0,1,false);  S=reshape(S(1,:,:),[1 input.cube.dsp.samples]);  psf=circshift(S,-23,2); edge=96;  psf=psf(1+edge:numel(psf)-edge);  psf=psf./sum(psf(:));
end

function psf=fcngetacrosspsf(input)
x=input.cube.pmt.stripx;  dx=x(2)-x(1);
s=input.cube.pmt.pdfacross; %4mm 1sigma charge cloud
psf = fcnpdfnorm((-5:5)*dx,0,s);   psf=psf/sum(psf(:));
end

function [x,y,z]=fcnfind2Dpeaks(X,threshold)
a = X<circshift(X,1,1);  a(1,:)=false;
b = X<circshift(X,1,2);  b(:,1)=false;

a = a & b & circshift(~a,1,1) & circshift(~b,1,2);
d=find(a);  d=d(X(d)>threshold);

[x,y]=fcnind2sub(size(a),d);  z=X(d);  n=numel(x); v=true(n,1);
[~,order]=sort(z,1,'descend'); %tallest maxima to lowest
rl = 11^2; %range limit 11 pixels
for i=order(:)'
    r = (x-x(i)).^2 + (y-y(i)).^2; %range squared
    j = r>0 & r<rl;  v(j)=false;  x(j)=nan;  y(j)=nan;
end
x=x(v);  y=y(v);  z=z(v);
end