% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = plotTSfiber()
clc; %close(findobj(0,'type','figure'));

file = uigetfile('/Users/glennjocher/Google Drive/MATLAB/neutrinos/nViewGUI/TSresults/*.mat'); fprintf('Loading ''%s''\n',file); %save uigetfilelast file
A=load(file);  MC = A.MC;    if ~iscell(MC.particleName); MC.particleName={MC.particleName}; end

%[C, XC, TC, EC, EFC, Cstr] = fcngetcandidates(A.input,MC.collectedPhotonCount,MC.xhat,MC.xtrue,MC.failure,'antineutrino');
n=numel(A.tsv);
for i=1
    j=min(i,numel(MC.particleName));
    plot1tsv(A,MC.particleName{j},MC.xtrue(:,:,i),MC.xhat(:,:,i),MC.vertex(:,:,i))
end

% if n>1
%     T=[MC.xtrue(:,28:29,1); MC.xtrue(:,28:29,3)];  X=[MC.xhat(:,28:29,1); MC.xhat(:,28:29,3)];
%     i = T(:,1) | T(:,2); figure; plotconfusion(T(i,:)',X(i,:)')
% end
end


function plot1tsv(A,particleName,T,X,vertex)
fs=12;
E = X-T;  EF = E./T; %error fraction
ec = 2; %edge cut (mm)
dt = (X(:,9)-X(:,4));
bi=~all(X==0,2) & ~any(isnan(X),2); %single scatter
switch particleName
    case {'FiberGamma','Gamma'}
        ci = bi & X(:,29)>.99 ...
            & X(:,14)>1 ... %> n recoils
            & all(abs(X(:,1:3))<(A.input.cube.Lr(3)-ec),2) ... %P1
            & all(abs(X(:,6:8))<(A.input.cube.Lr(3)-ec),2) ... %P2 
            & X(:,5)>.05  & X(:,5)<10 ... % E1
            & X(:,10)>.05 & X(:,10)<10 ... % E2
            & X(:,11) > 50 ... % P1 PE count
            & X(:,12) > 50 ... % P2 PE count
            & (X(:,5)>.2 | X(:,10)>.2) ...
            & fcnrange(X(:,1:3),X(:,6:8)) > 10 ... %dx (mm)
            & dt>.01 & dt<10 ... %dt (ns)
            & X(:,13)>1 & X(:,13)<179 ... %deg
            & X(:,15)<900 & X(:,15)>40 ... %V1 (mm/ns)
            & X(:,17)<10 & X(:,17)>0; %TRUE E& X(:,28)>.95
    case {'FiberNeutron','Neutron'}
        ci = bi ...            & X(:,28)>.99 ...
            & X(:,14)>1 ... %> number recoils
            & all(abs(X(:,1:3))<(A.input.cube.Lr(3)-ec),2) ... %P1 wall cut (mm)
            & all(abs(X(:,6:8))<(A.input.cube.Lr(3)-ec),2) ... %P2 wall cut (mm)
            & X(:,11)./X(:,12)>.1 ... % PE ratio
            & X(:,5)>.01 & X(:,5)<10 ... % E1 (MeV)
            & X(:,10)>.01 & X(:,10)<10 ... % E2 (MeV)
            & X(:,11) > 10 ... % P1 PE count (was 50)
            & X(:,12) > 10 ... % P2 PE count (was 50)
            & fcnrange(X(:,1:3),X(:,6:8)) > 10 ... %dx (mm)
            & dt>.1 & dt<20 ... %dt (ns) (was 1-20)
            & X(:,15)<40 & X(:,15)>5 ... %V1 (mm/ns)
            & X(:,17)<10 & X(:,17)>0; %TRUE E
end

fprintf('%.3f candidate efficiency\n',mean(ci));  %candidates
sb = 'true'; %----------------------------
switch sb
    case 'true';        e = T;
    case 'estimated';   e = X;
    case 'error';       e = E;
end
c1=fcndefaultcolors(1);  c2=fcndefaultcolors(2);  c3=fcndefaultcolors(3);


%POSITION RESOLUTION ------------------------------------------------------
[ha,hf]=fig(6,2,'54x18cm');  v=linspace(-1,1,100)*A.input.cube.Lr(3);
sca; fcnhistc(E(ci,3+0),v,c1,'bar','Z');  fcnhistc(E(ci,(2)+0),v,c2,'bar','Y');  fcnhistc(E(ci,(1)+0),v,c3,'bar','X');  xyzlabel(sprintf('%.3f efficiency {\\bf%s} P_1 (mm)',mean(ci),particleName)); alpha(.6)
sca; fcnhistc(E(ci,3+5),v,c1,'bar','Z');  fcnhistc(E(ci,(2)+5),v,c2,'bar','Y');  fcnhistc(E(ci,(1)+5),v,c3,'bar','X');  xyzlabel('P_2 (mm)'); alpha(.6)
sca; fcnhistc(E(ci,4+0),50,c1,'bar','T');  xyzlabel('T_1  (ns)');
sca; fcnhistc(E(ci,4+5),50,c1,'bar','T');  xyzlabel('T_2 (ns)');
sca; fcnhistc(E(ci,5+0),50,c1,'bar','E');  xyzlabel('E_1 (MeV)');
sca; fcnhistc(E(ci,5+5),50,c1,'bar','E');  xyzlabel('E_2 (MeV)');
sca; a=X(ci,14); fcnhistc(a,50);  xyzlabel('Recoils Observed')
sca; a=E(ci,17)./T(ci,17); [y,x]=fcnhistc(a,50,fcndefaultcolors(1));  xyzlabel('E_0 Resolution (Fraction)');  f=histNfit(x,y,'b');
sca; a=fcnsigmarejection(X(ci,15),3,6); fcnhistc(a,50);  xyzlabel('P_1-P_2 Velocity (mm/ns)')
sca; a=cosd(X(ci,16)); fcnhistc(a,50); xyzlabel(sprintf('Angle Error (cos\\theta), %.1f^\\circ 1\\sigma',acosd(1-fcnstdnonlin(a,.683,1))))
sca; [y,x]=fcnhistc(E(ci,13),50,c2,'bar'); histNfit(x,y,'k'); xyzlabel('Cone Error (\theta)');
sca; [y,x]=fcnhistc(X(ci,16),50,c2,'bar'); histNfit(x,y,'k'); xyzlabel('Angle Error (\theta)');
for i=1:numel(ha); axis(ha(i),'tight'); legend(ha(i),'show'); set(ha(i),'YTick',[],'YColor',[1 1 1]); grid(ha(i),'off'); end; hf.Position(1)=2.5; fcnfontsize(fs);  


try
    ha=fig(2,2,'19cm');  s='3L';
    x=vertex(:,5); %E
    %x=fcnel(vertex(:,6:8)); %el (rad)
    %x=vertex(:,3); %z (mm)
    sca(1); y=E(ci,17)./T(ci,17);              [mu,sigma,~,p]=movingMean(x(ci),y,100,0);            plot(p,sigma,'-','Display',sprintf('%s \\mu=%.3g',s,mean(sigma)));  xyzlabel('','Energy Resolution (fraction)')
    sca(2); y=X(ci,16);                        [mu,sigma,~,p]=movingMean(x(ci),y,100,0);            plot(p,sigma,'-','Display',sprintf('%s \\mu=%.3g',s,mean(sigma)));  xyzlabel('','Angle Resolution (deg)')
    sca(3); y=ci;                              [mu,sigma,~,p]=movingMean(x,y,100,0);                plot(p,mu,'-','Display',sprintf('%s \\mu=%.3g',s,mean(mu)));        xyzlabel('','Candidate Efficiency')
    %sca(4); y=1-abs(T(bi,28)-round(X(bi,28))); [mu,sigma,~,p]=movingMean(x(bi),y,100,0); plot(p,mu,'-','Display',sprintf('%s \\mu=%.3g',s,mean(mu)));        xyzlabel('','Discrimination')
    xyzlabel(ha,[particleName ' E (MeV)']); fcnlinewidth(2); fcntight('y0'); fcntight('x'); for i=1:4; legend(ha(i),'show'); end
end

% ''
% fig(2,2,'19cm');  x=vertex(:,5);
% sca; y=E(ci,17)./T(ci,17);          plotscatter(x(ci),y,'error (MeV)');
% sca; y=X(ci,16);                    plotscatter(x(ci),y,'error (deg)');
% sca; y=ci;                          plotscatter(x,y,'candidate efficiency');
% sca; y=1-abs(T(ci,28)-round(X(ci,28)));    plotscatter(x(ci),y,'discrimination');  xyzlabel('E (MeV)'); fcnlinewidth(2)

% %ANGULAR RESOLUTION -------------------------------------------------------
a = vertex(ci,6:8);
b = X(ci,6:8)-X(ci,1:3);
% % [SNR, ct] = fcnangles2SNR(b,a);
% % %fcnSNR2angle1sigma(.09,2500); %CHOOZ
b_a = rotateB2Wc(double(fcnVEC2RPY(b)),double(a)); %a expressed in b frame
% ct=cosd(X(ci,16)); %cos( abs(acos(ct)-d2r*X(ci,16)) );
%
% 
%PLOT SKYMAP --------------------------------------------------------------
a=cosd(X(ci,16));  nb=100;  x=linspace(-1,1,nb);  y=fcnhistc(a,x);
dx = 2/nb; ynormalized = y/sum(y)/dx/(2*pi);  F=fit(x(:),ynormalized(:),'spline'); %integrate to 1/(2pi)
%fig; bar(x,ynormalized); plot(x,F(x))
for i=[1 3 10 100]
    fig(1,2,2);  plotSkyMap(F,b_a,real(X(ci,13)),i,particleName);
end



% [ha,hf]=fig(5,2,'45x18cm');  ct=E(ci,16);
% sca; a = (X(ci,10)-X(ci,5));                            plotscatter(a,ct,'cos(\theta)'); xyzlabel('E_1 (MeV)');
% sca; a = X(ci,5)./X(ci,10);                             plotscatter(a,ct,'cos(\theta)'); xyzlabel('E_2 (MeV)'); 
% sca; a = fcnelaz(vertex(ci,6:8))*r2d; a=a(:,1);         plotscatter(a,ct,'cos(\theta)'); xyzlabel('el (deg)');
% sca; a = X(ci,3);                                       plotscatter(a,ct,'cos(\theta)'); xyzlabel('P_1Z (mm)'); 
% sca; a = rangec(double(X(ci,6:8)-X(ci,1:3)));           plotscatter(a,ct,'cos(\theta)'); xyzlabel('dx (mm)');
% sca; a=X(ci,9)-X(ci,4);                                 plotscatter(a,ct,'cos(\theta)'); xyzlabel('dt (ns)');
% sca; a=T(ci,18);                                        plotscatter(a,ct,'cos(\theta)'); xyzlabel('True E_0 (MeV)');
% sca; a=X(ci,15);                                        plotscatter(a,ct,'cos(\theta)'); xyzlabel('V (mm/ns)');
% sca; plotscatter(vertex(:,5),ci,'candidate'); xyzlabel('E_0 (MeV)');
% sca; plotscatter(X(ci,13),E(ci,17),'cos(\theta)'); xyzlabel('V (mm/ns)');
% fcntight;  fcnfontsize(fs); hf.Position(1)=22.5; %fcntight(ha,'csigma'); % set(ha,'ylim',[-1 1]);

%find(X(:,14)==3 & E(:,17)>3 & ci) %HIGH ERROR EVENTS
 

% [ha,hf]=fig(5,2,'45x18cm'); L=A.input.cube.Lr;
% for i=1:2
%     j=(i-1)*5;
%     %sca(ha(0+i)); a=T(ci,1+j); b=X(ci,1+j);     plotscatter(a,b); xyzlabel(sprintf('X_%g (mm)',i),'fit (mm)');
%     %sca(ha(2+i)); a=T(ci,2+j); b=X(ci,2+j);     plotscatter(a,b); xyzlabel(sprintf('Y_%g (mm)',i),'fit (mm)');  
%     %sca(ha(4+i)); a=T(ci,3+j); b=X(ci,3+j);     plotscatter(a,b); xyzlabel(sprintf('Z_%g (mm)',i),'fit (mm)');   plot([-100 100],[-100 100],'w-','Linewidth',1)
%     sca(ha(0+i)); a=T(ci,1+j); b=X(ci,1+j);     xb=linspace(-L(1),L(1),100); histogram2(a,b,xb,xb,'DisplayStyle','Tile');  fcntight(gca,'csigma');  xyzlabel(sprintf('X_%g (mm)',i),'fit (mm)');
%     sca(ha(2+i)); a=T(ci,2+j); b=X(ci,2+j);     xb=linspace(-L(2),L(2),100); histogram2(a,b,xb,xb,'DisplayStyle','Tile');  fcntight(gca,'csigma');  xyzlabel(sprintf('Y_%g (mm)',i),'fit (mm)');
%     sca(ha(4+i)); a=T(ci,3+j); b=X(ci,3+j);     xb=linspace(-L(3),L(3),100); histogram2(a,b,xb,xb,'DisplayStyle','Tile');  fcntight(gca,'csigma');  xyzlabel(sprintf('Z_%g (mm)',i),'fit (mm)');
%     
%     sca(ha(6+i)); a=T(ci,4+j); b=X(ci,4+j);     plotscatter(a,b); xyzlabel(sprintf('T_%g (ns)',i),'fit (ns)');   plot([0 max(a)],[0 max(a)],'w-','Linewidth',1)
%     sca(ha(8+i)); a=T(ci,5+j); b=X(ci,5+j);     plotscatter(a,b); xyzlabel(sprintf('E_%g (MeV)',i),'fit (MeV)'); plot([0 max(a)],[0 max(a)],'w-','Linewidth',1)
% end
% fcntight;   for k=1:5;  ha(k*2).XLim=ha(k*2-1).XLim;  ha(k*2).YLim=ha(k*2-1).YLim;  end;  fcnfontsize(fs);  hf.Position(1)=42.5;
% 
% [~,hf]=fig(5,2,'45x18cm');
% sca; a=T(ci,13); b=X(ci,13);	plotscatter(a,b); xyzlabel('angle (deg)','fit');
% sca; a=T(ci,15); b=X(ci,15);	plotscatter(a,b); xyzlabel('V (mm/ns)','fit');
% sca; a=T(ci,17); b=X(ci,17);	plotscatter(a,b); xyzlabel('E_0 (MeV)','E_0 corrected (MeV)');
% sca; a=T(ci,18); b=X(ci,18);	plotscatter(a,b); xyzlabel('E_0 (MeV)','E_0 uncorrected (MeV)');
% sca; a=X(ci,14); b=E(ci,17);	plotscatter(a,b); xyzlabel('N Recoils','E_0 corrected (MeV)');
% sca; a=X(ci,14); b=E(ci,18);	plotscatter(a,b); xyzlabel('N Recoils','E_0 uncorrected (MeV)');
% fcntight;  fcnfontsize(fs);   hf.Position(1)=62.5;
end

function c2() %prompt color
h=findobj(gca,'type','line');  h.Color = fcndefaultcolors(2);
end
function c3() %delayed color
h=findobj(gca,'type','line');  h.Color = fcndefaultcolors(3);
end

function plotscatter(x,a,yl)
if nargin<3; yl=''; end

h=histogram2(x,a,100,'DisplayStyle','Tile');  movingMean(x,a,100,1); 
xyzlabel('E_\nu (MeV)',yl); %legend('Location','SouthEast')
%h.DisplayName=sprintf('R=%.2f',corr(x,a));
colorbar off
end

function plothist(x,a,dn,yl)
if nargin==3; yl=''; end
[hy,hx]=fcnhistc(a,50,[],'bar',dn);  xlabel(yl); histNfit(hx,hy,'b'); legend show; 
end
