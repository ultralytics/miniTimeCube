% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnfindplutonium2()
clc; close all
%http://www-pub.iaea.org/MTCD/publications/PDF/P1433_CD/datasets/abstracts/sm_en-20.html
[d, s, DEM] = fcngetpositions();



load Goldhagenneutron.mat
yb = y;
xb = x;

load Plutoniumneutron.mat
xs = e;  dxs = e(2)-e(1);
ys = pdfs; %y signal

B = load('MC Neutrons-2L Goldhagen Background PLANACON XP85012-A1 80pc 1536pix 20qe 13re 2cp 1%EJ-254 10k Neutron.mat');
B.z = fcngetcandidates(B.MC);
B.epdfx = xs;
B.epdfy = fcnanalyticspectrum(xs+eps);
B.epdfy = B.epdfy/sum(B.epdfy)/dxs;

S = load('MC Neutrons-2L Pu Spectrum PLANACON XP85012-A1 80pc 1536pix 20qe 13re 2cp 1%EJ-254 10k Neutron.mat');
S.z = fcngetcandidates(S.MC);
S.epdfx = xs;
S.epdfy = ys;
S.epdfy = S.epdfy/sum(S.epdfy)/dxs;




h=fig(2,2,2);
sca(h(1))
plot(S.epdfx,S.epdfy,'b'); hold on; plot(B.epdfx,B.epdfy,'r'); 
legend('Pu Spectrum','Background Spectrum'); xlabel('Neutron Energy (MeV)'); ylabel('pdf'); grid on


sca(h(2))
x1 = linspace(0,10,30);
y1 = hist(B.z.ehat,x1); bar(x1,y1,1); axis tight
sca(h(3))
y1 = hist(S.z.ehat,x1); bar(x1,y1,1); axis tight



t = 3600;
fluxb = .0134*169*3*.065;
nd = numel(d);
for i=1:nd
    d(i).flux(2) = fluxb;
    d(i).nmean = d(i).flux*t;
    d(i).nrandom = [fcnRandomPoisson(d(i).nmean(1)) fcnRandomPoisson(d(i).nmean(2))];
    
    nj = d(i).nrandom(1);  vnj = repmat(1:S.z.n,[1, ceil(nj/S.z.n)+10]);  j = vnj(1:nj);
    nk = d(i).nrandom(2);  vnk = repmat(1:B.z.n,[1, ceil(nk/B.z.n)+10]);  k = vnk(1:nk);
    
    d(i).z.ehat = [S.z.ehat(j); B.z.ehat(k)];
    d(i).z.elhat = [S.z.elhat(j); B.z.elhat(k)];
    d(i).z.azhat = [S.z.azhat(j); B.z.azhat(k)];
    d(i).z.vechat = [S.z.vechat(j,:); B.z.vechat(k,:)];

    
    d(i).est.rmat = sqrt((DEM.estx - d(i).pos(1)).^2 + (DEM.esty - d(i).pos(2)).^2);

    m = 1; %kg of pu
    n = 6.1E4*m; % (n/s) Weapons Grade Plutonium WGP = 6% Pu240, 94% Pu239
    area = 13^2; %detector area (cm^2)
    flux = n./(4*pi*(d(i).est.rmat*100).^2)*area; %flux (n/m^2/s)
    efficiency = .12;
    d(i).est.fluxperkg = flux*efficiency + fluxb;
end




























%PU FLUX
m = 1; %kg of pu
n = 6.1E4*m; % (n/s) Weapons Grade Plutonium WGP = 6% Pu240, 94% Pu239
r = 7; % (m) range
area = 13^2; %detector area (cm^2)
f = n/(4*pi*(r*100)^2); %flux (n/m^2/s)
ys = ys*f;



e1 = .1; %(MeV) energy floor
e2 = 10.2; %(MeV) energy ceiling

ybs = yb * area; %y background scaled
yss = ys * area; %y signal scaled

ib = xb>e1 & xb<e2;
is = xs>e1 & xs<e2;
xa = xb(ib);
yba = fcnanalyticspectrum(xa);

h=fig(2,1,2,1);
sca(h(1))
plot(xb(ib),yb(ib)./xb(ib),'b.-'); hold on
plot(xs(is),ys(is),'r.-')
plot(xa,yba,'g.-'); axis tight
fa = trapz(xa,yba*area);
fb = trapz(xb(ib),ybs(ib)) %flux background (n/s)
fs = trapz(xs(is),yss(is)); %flux signal (n/s)
xlabel('MeV'); ylabel('n/cm^2/s'); axis tight; set(gca,'xlim',[0 e2]);
legend(sprintf('%.3f/cm^2/s   Sea Level Neutrons (Goldhagen 2004)',  fb/area), ...
    sprintf('%.3f/cm^2/s   1kg Plutonium Neutrons at %.1fm (Oshawa 2000)',  fs/area,r), ...
    sprintf('%.3f/cm^2/s   f(E)/E = 1e-6exp(2.1ln(E)-0.35ln(E)^2) + 1e-3exp(-0.67ln(E)-0.41ln(E)^2)',  fa/area));

sca(h(2))
semilogx(xb,yb,'b.-'); hold on
plot(xs(is),ys(is),'r.-')
plot(xb,fcnanalyticspectrum(xb).*xb,'g.-');
xlabel('MeV'); ylabel('n/cm^2/s'); axis tight; title('Semilogx Version. Compare to Goldhagen Figures 6 and 7.')

%integral(@(x)fcnanalyticspectrum(x),e1,e2)

[~,i] = fcnunique(xb); xb=xb(i); yb=yb(i); xb=log10(xb);

xb2 = linspace(min(xb),max(xb),1000);
yb2 = interp1(xb,yb,xb2,'linear');

xx = fcnrandcdf(cumsum(yb2),xb2,10000,'linear');  %xx=10.^xx;
pd = fitdist(xx,'Kernel','Bandwidth',.1);
%plot(10.^xb2,pdf(pd,xb2)/175,'c')

figure; hist(xx,50)
figure; plot(xb2, pdf(pd,xb2));

%xx = fcnrandcdf(cumsum(pdf(pd,xb2(1:793))),xb2(1:793),10000,'linear');  xx=10.^xx;  save xx.mat xx
%fa = trapz(xb,yb)
end

function [d, s, DEM] = fcngetpositions()
%GET GE IMAGE -------------------------------------------------------------
%cam.ssge = [8 2560-720-77 1280 720];  %1080x720; 27" 27" dual monitor hack
%he = actxserver('googleearth.ApplicationGE'); %Create a COM server running Google Earth
%DEM = fcngetGEDEM(he,70);  fcnFinishStreamProgress(he);  x=getscreen(cam.ssge);  DEM.cdata = x.cdata;  save DEMHamburg.mat DEM
%MOL Cosmos, Panama

load DEMHamburg
fig(1,1,3.8,2.5);
sd = size(DEM.cdata);

pix2m = (abs(diff(DEM.nedx(1,[1 end])))/sd(1) + abs(diff(DEM.nedy([1 end],1)))/sd(2))/2*1E3;
ax = linspace(1,sd(1),sd(1))*pix2m; %up-down
ay = linspace(1,sd(2),sd(2))*pix2m; %left-right
imshow(DEM.cdata,'Xdata',ay,'Ydata',ax); hold on; axis on; box off
title('Place Plutonium Source'); xlabel('meters'); ylabel('meters')

n = size(DEM.lat,1);
[DEM.y, DEM.x] = ndgrid(linspace(ax(end),ax(1),n),linspace(ay(1),ay(end),n));

for i=1:30
    [xi,yi,b] = ginput(1);
    if b~=1;break; end
    switch i
        case 1
            plot(xi,yi,'o','MarkerSize',10,'MarkerEdgeColor',[1 1 1]*.9,'MarkerFaceColor','r');
            text(xi,yi,sprintf('  1kg Plutonium'),'color',[1 1 1]*.9,'fontsize',8);
        otherwise
            plot(xi,yi,'o','MarkerSize',6,'MarkerEdgeColor',[1 1 1]*.9,'MarkerFaceColor','g');
            ht(i-1) = text(xi,yi,sprintf('  D%.0f',i-1),'color',[1 1 1]*.9,'fontsize',8);
    end
    x(i)=xi; y(i)=yi;
    title(sprintf('Place Detector Number %.0f, or right click to finish',i))
end
nd = size(x,2)-1;

lla(:,1) = interp2(DEM.x,DEM.y,DEM.lat',x,y);  lla(:,2) = interp2(DEM.x,DEM.y,DEM.lng',x,y);
lla(:,3) = DEM.F(lla(:,2),lla(:,1)); %m above ellipsoid
ecef = lla2ecef(lla);

s.pos = [x(1) y(1)];
s.poslla = lla(1,:);
s.posecef = ecef(1,:);
for i=1:nd
    j = i+1;
    d(i).pos = [x(j) y(j)]; %#ok<*AGROW>
    d(i).poslla = lla(j,:);
    d(i).posecef = ecef(j,:);
    
    r = fcnrange(d(i).posecef,s.posecef)*1E3; %m
    m = 1; %kg of pu
    n = 6.1E4*m; % (n/s) Weapons Grade Plutonium WGP = 6% Pu240, 94% Pu239
    area = 13^2; %detector area (cm^2)
    flux = n/(4*pi*(r*100)^2)*area; %flux (n/m^2/s)
    efficiency = .12;
    
    d(i).r = r; %m
    d(i).flux = flux * efficiency; %n/cm^2/s
    
    set(ht(i),'string',sprintf('  D%.0f\n  %.0fm\n  %.3fn/s',i,r,flux))
end

%Estimation Map
[DEM.esty, DEM.estx] = ndgrid(linspace(ax(end),ax(1),200),linspace(ay(1),ay(end),300));
end

function z = fcngetcandidates(MC)
    i = 1;
    a = zeros(size(MC.xhat(:,1,1),1),13);
    a(:,1) = MC.xhat(:,1,i)~=0;
    a(:,2) = MC.xhat(:,6,i)>10; %P1 > 10 photons
    a(:,3) = MC.xhat(:,12,i)>5; %P2 > 2  photons
    a(:,4) = MC.xhat(:,10,i) - MC.xhat(:,4,i) > 1; %dt > 1ns
    a(:,5) = fcnrange(MC.xhat(:,1:3,i) - MC.xhat(:,7:9,i)) > 10; %dx > 10mm
    a(:,6) = MC.xhat(:,5,i)>0.200; %.3 %first bounce > 100keV
    a(:,7) = MC.xhat(:,6,i)./MC.xhat(:,12,i) > .2; %dE1/dE2 must be greater than 0.20
    a(:,8) = MC.xhat(:,5,i)./MC.xhat(:,18,i) < .9; %dE1./E0 < 0.9
    a(:,9) = MC.xhat(:,5,i)./MC.xhat(:,18,i) > .1; %dE1./E0 > 0.1
    a(:,10) = MC.xhat(:,5,i)<3; %.3 %first bounce > 100keV
    %a(:,10) = fcninsidedetector(MC.xhat(:,1:3,i), input, 10); %P1 wall cut
    %a(:,11) = fcninsidedetector(MC.xhat(:,7:9,i), input, 10); %P2 wall cut    
    
    %a(:,12) = MC.collectedPhotonCount(:,2,i)>2; %delayed bookend
    %a(:,13) = all(MC.protonprotonflag(:,:,i),2); %proton-proton bounces only!
    j = all(a(:,[1:10]),2);    sum(j)/numel(j)
    nj = sum(j);
    
    
    z.n = nj;
    z.ehat = MC.xhat(j,18,i); 
    p1 = MC.xhat(j,1:3,i);
    p2 = MC.xhat(j,7:9,i);
    z.vechat = fcnvec2uvec(p2-p1);
    sc = fcnCC2SC(z.vechat);
    z.elhat = sc(:,2);
    z.azhat = sc(:,3);
    z.anglehat = MC.xhat(j,13,i);



    %MC.xhat(j,19,i) = input.neutron.dEfit(MC.xhat(j,19,i));
    %p1 = MC.xhat(j,1:3,i);
    %p2 = MC.xhat(j,7:9,i);
    %pn = MC.xhat(j,15:17,i);
    %anglehat = MC.xhat(j,13,i);
    %E0hat = MC.xhat(j,18,i);
end


function f = fcnanalyticspectrum(E)
% %http://phits.jaea.go.jp/expacs/data/Sato-RR166-p544-2006.pdf
% c1 = .229;
% c2 = 2.31; 
% c3 = .721;
% c4 = .0516;
% c5 = 126;
% c6 = 2.17;
% c7 = .00108;
% c8 = 3.33E-12;
% c9 = 1.62;
% c10 = 9.59E-8;
% c11 = 1.48;
% c12 = 299;
% f = c1*(E/c2).^c3.*exp(-E/c2)+c4*exp((-(log10(E) - log10(c5)).^2)./(2*log10(c6).^2)) ...
%     + c7*log10(E/c8).*(1 + tanh(c9*log10(E/c10))).* ...
%     (1 - tanh(c11*log10(E/c12)));
% f = f/100;

%Goldhagen analytic EQN6
B1 = .3500;
B2 = .4106;
y1 = 2.1451;
y2 = -.6670;
c1 = 1.006E-6;
c2 = 1.011E-3;

f =  c1*exp(-B1*log(E).^2 + y1*log(E)) ...
   + c2*exp(-B2*log(E).^2 + y2*log(E));

end