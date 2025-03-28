% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnPlotWavelengthDepDist(input)
%Plots wavelength dependent distributions
ha = fig(2,3);
%ha = fig(3,1,'27x18cm');
%ha = fig(1,3,'10x30cm');

F = input.Material(1);
X = input.wl;
size = 40; 
xlim = [201 849];
clight=299.792458; % (mm/ns)

%\lambda Probability Distribution
sca
wl = input.wl;  v1 = 30:fcnindex1c(wl,X(end));  wl=wl(v1);
scatter(X, F.X(:,5), size, X, 'filled') 

ylabel('pdf')
title([input.Material(1).name ' Scintillation Spectrum'])
colormap(input.colormap)
h = colorbar('location','East');  h.Limits = xlim;  set(h,'YTick',300:100:700,'FontSize',9);
fcncolorbar(1,'nm'); h.Box='off';

%cherenkov
sca
wlpdf = cherenkov_magnitude(wl, F.X(v1,1), 1); wlpdf=wlpdf/sum(wlpdf);
scatter(wl, wlpdf, size, wl, 'filled')
%scatter(in.x, in.wlpdf, size, in.x, 'filled')
ylabel('pdf')
title('Cherenkov Spectrum')

sca
scatter(X, F.X(:,2), size, X, 'filled')
title('Attenuation Length (SG BCF-98)')
ylabel('Attenuation Length (m)')

sca
scatter(X, F.X(:,1), size, X, 'filled')
title('n (Sellmeier: B1=1.44, C1=82.8nm^2)')
ylabel('Index of Refraction');
%ylabel('photon speed (mm/ns)'); ha(4).YLim = [187 192];

sca
scatter(X, F.X(:,4), size, X, 'filled')
title('Re-Emission Efficiency (Oleg)')
ylabel('Re-Emission Efficiency')

sca
scatter(X, input.cube.QE, size, X, 'filled')
title([input.cube.prettyname ' QE'])
ylabel('Quantum Efficiency')

fcntight(ha,'xy')
xyzlabel(ha,'\lambda (nm)')

ha(2).YLim(2) = 6E-3;
%ha(4).YLim(2) = 1.75;
set(ha,'xlim',xlim);
fcnfontsize(14)

%return

%PHOTON PHASE AND GROUP VELOCITIES!
ha = fig(2,2,'27x27cm');
sca
scatter(X, F.X(:,5), size, X, 'filled') 
ylabel('pdf')
title([input.Material(1).name ' Scintillation Spectrum'])
colormap(input.colormap)
h = colorbar('location','East');  h.Limits = xlim;  %set(h,'YTick',200:100:700,'FontSize',9);
fcncolorbar(1,'nm'); h.Box='off';

sca
scatter(X, F.X(:,1), size, X, 'filled')
title('n (Sellmeier: B1=1.44, C1=82.8nm^2)')
ylabel('Index of Refraction');

sca
speed=clight./F.X(:,1);  scatter(X, speed, size, X, 'filled','Display','Phase Velocity')
dlambda=1;  dn=[0; diff(F.X(:,1))];  speed=clight./(F.X(:,1) - X.*dn./dlambda);  scatter(X, speed, size, X, 'filled','Display','Group Velocity')
ylabel('photon speed (mm/ns)');
xyzlabel(ha,'\lambda (nm)')

sca
wl = randcdfc(input.Material(1).R.scintWL,1E6);
G = griddedInterpolant(X,speed);  speed=G(wl);
fhistogram(fcnsigmarejection(speed),30); xlabel('group speed (mm/ns) = c/(n-\lambda*dn/d\lambda)')

fcntight; ha(2).YLim = [1.55 1.7];  ha(3).YLim = [170 192];  set(ha(1:3),'XLim',[200 800])
a=1000./speed;

