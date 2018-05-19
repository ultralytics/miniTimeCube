clc; clear all; close all
load Goldhagenneutron.mat

% MWth=20;
% y=y./x;
% plot(x,y,'r','linewidth',2); axis tight; grid on
% loglog(x,y,'r','linewidth',2); axis tight; grid on
% %trapz(x,y.*x)


dx = x(end)-x(1);
y1 = y./x;
fig; semilogy(x,y1); hold on

load Plutoniumneutron.mat
pdfs = pdfs/100;
plot(e,pdfs,'r')
xlabel('MeV'); ylabel('n/cm^2/s');
fcnfontsize(8)
str{1}='Background Neutrons (Goldhagen 2004)';
str{2}='Plutonium Neutrons (Oshawa 2000)';
axis tight
set(gca,'xlim',[0 50],'ylim',[1E-5 max(pdfs)]);

%BACKGROUND FLUX
sum(y1(x<.001)) %n/cm^2/s between 1 and 11 MeV!

%PU FLUX
n = 6.4E4; % (n/s) Weapons Grade Plutonium WGP = 6% Pu240, 94% Pu239
r = 25; % (m) range
f = 6.1E4/(4*pi*r^2); %flux (n/m^2/s)

%REACTOR NEUTRONS ---------------------------------------------------------
clear pdf
x = linspace(0,20,2000);
k = 1.5; 
theta = 0.64;
y = pdf('Gamma',x,k,theta);
plot(x,y,'g');
p = sum(y(1:199))/sum(y)*100; %percent under 2MeV
str{3} = sprintf('Reactor Neutrons: gamma pdf, k=%3g, \\theta=%3g (%.1f%%<2MeV)',k,theta,p);

legend(str)
set(gca,'yscale','log')
set(gca,'xlim',[0 20])

% %tricolor shawn plot
% fig; i=1:48; area(log10(x(i)),y(i),0,'facecolor',[1 0 0],'edgecolor','none'); i=48:90; area(log10(x(i)),y(i),0,'facecolor',[0 1 0],'edgecolor','none'); i=90:129; area(log10(x(i)),y(i),0,'facecolor',[0 0 1],'edgecolor','none')
% xyzlabel('log neutron energy (log MeV)',sprintf('Energy x Differential Flux\n E dPhi//dE (cm^-2s^-1)')); grid on
% fcnfontsize(12)