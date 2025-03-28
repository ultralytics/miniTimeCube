% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = plotTScontinuous()
clc; close(findobj(0,'type','figure'))

file = uigetfile();
A=load(file);  MC = A.MC;

[C, XC, TC, EC, EFC, Cstr] = fcngetcandidates(A.input,MC.collectedPhotonCount,MC.xhat,MC.xtrue,MC.failure,'antineutrino');
T = MC.xtrue;  X = MC.xhat;  if min(T(:,4))>1.78; T(:,4)=T(:,4)-1.807; end
E = X-T;  EF = E./T; %error fraction
PE = double(MC.collectedPhotonCount);  ci = C(:,end);
x = T(:,4)+1.8 + randn(size(T(:,4)))*.01; %x axis

% ha=fig(1,1,'18x18cm'); [xs,j]=sort(x);
% for k=1:numel(Cstr)
%     a = fcnsmooth(C(j,k),399);  plot(xs,a,'.','displayname',Cstr{k},'tag','a')
% end; fcnmarkersize(20); fcntight; legend show; xyzlabel(ha,'E_\nu (MeV)','candidates'); fcnfontsize(16); %set(gca,'xlim',[2 9.5]);
% fcnplotsigmabounds(.683,30); delete(findobj(gcf,'tag','a')); fcnfontsize(16);


sa = 'candidates'; %----------------------------
switch sa
    case 'all';         i=true(size(x));
    case 'candidates';  i=ci;
end
xi=x(i);

sb = 'true'; %----------------------------
switch sb
    case 'true';        e = T;
    case 'estimated';   e = X;
    case 'error';       e = E;
end

sc = 'delayed'; %----------------------------
switch sc
    case 'prompt';      m=0; pei=1;
    case 'delayed';     m=5; pei=2;
end
dn = [sb ' ' sc ' ' sa];

% %GENERAL ------------------------------------------------------------------
% ha=fig(4,1,1.5);
% sca; a = rangec(e(i,(1:3)+m));      plotscatter(xi,a,dn,'P (mm)')
% sca; a = e(i,4+m);                  plotscatter(xi,a,dn,'E (MeV)')
% sca; a = e(i,5+m);                  plotscatter(xi,a,dn,'T (ns)')
% sca; a = PE(i,pei);                 plotscatter(xi,a,dn,'PE')
% fcntight(ha,'xy'); fcnfontsize(16);  fcnmarkersize(5);  fcnlinewidth(2); fcntight('y sigma')


% %POSITRON -----------------------------------------------------------------
% ha=fig(2,2,'18x18cm');  %xi=X(i,1);  %VS NUEBAR ENERGY
% sca;                                                sca; plotscatter(xi,rangec(E(i,1:3)),'point fit','P (mm)'); 
% sca; plotscatter(xi,X(i,26),'e^+ fit','cos\theta'); sca; plotscatter(xi,rangec(E(i,23:25)),'e^+ fit','P (mm)'); 
% fcnfontsize(16);  fcnlinewidth(2); fcntight('xjoint'); fcntight(ha([2 4]),'ysigma joint')

%XY POSITION SCATTER 
fig(2,2);
sca; hist2(T(ci,1),X(ci,1),50,'pcolor');  xyzlabel('X_{true} (mm)','X_{est} (mm)','','Prompt')
sca; hist2(T(ci,6),X(ci,6),50,'pcolor');  xyzlabel('X_{true} (mm)','X_{est} (mm)','','Delayed')
sca; hist2(X(ci,1),X(ci,2),50,'pcolor');  xyzlabel('X_{est} (mm)','Y_{est} (mm)','','Prompt')
sca; hist2(X(ci,6),X(ci,7),50,'pcolor');  xyzlabel('X_{est} (mm)','Y_{est} (mm)','','Delayed');  fcntight
fcnfontsize(16);  fcntight('x y sigma');  

%ENERGY RESOLUTION --------------------------------------------------------
fig; a = EF(ci,4);  plot(x(ci),a,'.'); [~,s,xb]=fcnplotsigmabounds(.683,30); delete(gcf);
ha=fig(2,1,'18x18cm'); %VS NUEBAR ENERGY
sca; plotscatter(x(ci),a,'prompt energy error'); ylabel('(Est-True)/True (fraction)')
sca; plot(xb,s,'.-'); ylabel('E_\sigma (fraction)');  xyzlabel(ha,'E_\nu (MeV)');  
fcnfontsize(16);  fcnmarkersize(ha(2),25);  fcnlinewidth(2); fcntight('ysigma'); fcntight(ha(2),'y0'); fcntight('xjoint')

fig; plot(X(ci,1),a,'.'); [~,s,xb]=fcnplotsigmabounds(.683,30); delete(gcf);
ha=fig(2,1,'18x18cm'); %VS VERTEX LOCATION
sca; plotscatter(X(ci,1),a,'prompt energy error','E (MeV)'); ylabel('(Est-True)/True (fraction)')
sca; plot(xb,s,'.-'); ylabel('E_\sigma (fraction)');  xyzlabel(ha,'detector side (mm)');  
fcnfontsize(16);  fcnmarkersize(ha(2),25);  fcnlinewidth(2); fcntight('ysigma'); fcntight(ha(2),'y0'); ha(1).XLim = [-67 67]; ha(2).XLim = [-67 67];

%POSITION RESOLUTION ------------------------------------------------------
ha=fig(2,1,'18x18cm');
sca; a = rangec(E(ci,(1:3)+0));     plotscatter(x(ci),a,'prompt candidates','error (mm)');      %  fcnplotsigmabounds(.683,30);
sca; a = rangec(E(ci,(1:3)+5));     plotscatter(x(ci),a,'delayed candidates','error (mm)'); %c2; %  fcnplotsigmabounds(.683,30);
fcnfontsize(16);  fcntight('x y sigma'); ha(1).YLim(2)=40; ha(2).YLim(2)=40; 

%PE COUNTS ----------------------------------------------------------------
fig(2,1,'18x18cm'); %VS NUEBAR ENERGY
sca; a = PE(:,1);                   plotscatter(x,a,'all prompt','PE')
sca; a = PE(:,2);                   plotscatter(x,a,'all delayed','PE'); %c2
fcntight('xy sigma'); fcnfontsize(16); fcnlinewidth(2);

ha=fig(2,1,'18x18cm'); %VS VERTEX LOCATION
j = abs(T(:,6))<67;
sca; a = PE(:,1);                   plotscatter(T(:,1),a,'all prompt','PE');      %[~,s,xb]=fcnplotsigmabounds(.683,30);
sca; a = PE(j,2);                   plotscatter(T(j,6),a,'all delayed','PE'); %c2; %[~,s,xb]=fcnplotsigmabounds(.683,30);
fcntight('xy sigma'); fcnfontsize(16); fcnlinewidth(2); xyzlabel(ha,'detector side (mm)');  


%TIME AND DISTANCE WINDOW -------------------------------------------------
fig(2,1,'18x18cm');  c=fcndefaultcolors(1);
sca; a = (T(:,10)-T(:,5))/1E3;          fcnhistc(a,50,c,'bar',['all neutron dt ' sa]);  xyzlabel('neutron dt (\mus)');  
sca; a = rangec(T(:,6:8)-T(:,1:3));     fcnhistc(a,50,c,'bar',['all neutron dx ' sa]);  xyzlabel('neutron dx (mm)');  
fcntight('xy sigma'); fcnfontsize(16);  fcnmarkersize(25);  fcnlinewidth(2); set(findobj(gcf,'type','line'),'color',c);

%ANGULAR RESOLUTION ------------------------------------------------------
a = MC.vertex(ci,4:6);
b = X(ci,6:8)-X(ci,1:3);
weights = fcnspecreactoremission(x(ci)).*fcnspecdetection(x(ci));
[SNR, ct] = fcnangles2SNR(b,a,weights);
%fcnSNR2angle1sigma(.09,2500); %CHOOZ

ha=fig(3,2,'27x18cm');
sca; a = X(ci,4);                       plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('E_\nu (MeV)'); %fcnplotsigmabounds(.683,30);
sca; a = X(ci,9);                       plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('E_n (MeV)'); %c2; fcnplotsigmabounds(.683,30);
sca; a = X(ci,1);                       plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('P_\nu (mm)'); %fcnplotsigmabounds(.683,30); 
sca; a = X(ci,6);                       plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('P_n (mm)'); %c2; fcnplotsigmabounds(.683,30); 
sca; a = rangec(X(ci,6:8)-X(ci,1:3));   plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('dx (mm)'); %c3; fcnplotsigmabounds(.683,30);
sca; a = (X(ci,10)-X(ci,5))/1E3;        plotscatter(a,ct,sprintf('R=%.2f',corr(a,ct)),'cos(\theta)'); xyzlabel('dt (ns)'); %c3; fcnplotsigmabounds(.683,30);
fcntight; set(ha,'ylim',[-1 1]); set(ha(3:4),'xlim',[-67 67]); fcnfontsize(16)


%TIME AND DISTANCE WINDOW -------------------------------------------------
%fig; hist2(xi,ct,30,'pcolor');
end

function c2() %prompt color
h=findobj(gca,'type','line');  h.Color = fcndefaultcolors(2);
end
function c3() %delayed color
h=findobj(gca,'type','line');  h.Color = fcndefaultcolors(3);
end

function plotscatter(x,a,dn,yl)
if nargin==3; yl=''; end
%h=plot(x,a,'.','markersize',2);  
h=hist2(x,a,[100 50],'pcolor'); 
xyzlabel('E_\nu (MeV)',yl); legend show
h.DisplayName=dn;
colorbar off
end

function plothist(x,a,dn,yl)
if nargin==3; yl=''; end
[hy,hx]=fcnhistc(a,50,[],'bar',dn);  xlabel(yl); histNfit(hx,hy,'b'); legend show; 
end