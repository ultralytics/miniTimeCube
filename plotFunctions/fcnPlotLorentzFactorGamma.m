function fcnPlotLorentzFactorGamma(input, photons)
h=fig(1,2,1.5);
sca(h(1))
m = 0.510998910; %MeV
ke = linspace(0,3,1000); %MeV
vel = sqrt(ke.^2 + 2*m*ke)./(ke+m);

%PLOT SPEED VS KE
plot(ke,vel);
xlabel('Positron Kinetic Energy (MeV)')
ylabel('Positron Speed (c)')
title({'Positron Speed Using Lorentz Factor Gamma Equation','v = sqrt(ke^2 + 2m*ke) / (ke+m)'})

%PLOT MEAN SPEED OF PHOTONS IN THE MEDIUM
velPhotons = 1/mean(photons.refractiveIndex);
plot([min(ke) max(ke)],[1 1]*velPhotons,'r-','LineWidth',2)

%ADD LEGEND
legend('Positron Speed','Mean Photon Speed')

%PLOT CHERENKOV HALF ANGLE VS SPEED
sca(h(2))
vel = linspace(velPhotons,1,1000);
pitch = real(acosd(1./vel./input.Material(1).mu(1))); %cherenkov angle for a positron
plot(vel,pitch,'g-');
ylabel('Cherenkov Half-Angle Cone (deg)')
xlabel('Positron Speed (c)')
title({'Cherenkov Angle vs Positron Speed','theta = acos(1/(speed*refractiveIndex))'})
set(gca,'XDir','Reverse')

%PLOT MEAN SPEED OF PHOTONS IN THE MEDIUM
plot([1 1]*velPhotons,get(gca,'YLim'),'r-','LineWidth',2)
legend('Cherenkov Angle','Mean Photon Speed')

fcngrid(h,'on')
