% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function c = fcnke2c(ke,mass)
% kE in MeV
% mass in MeV
% m = 0.510998910; %e+ MeV
c = sqrt(ke.^2 + 2*mass*ke)./(ke+mass); %speed of particle (c fraction)


% ke = fcnc2ke(c,mass);
% mass = 939.565; %(MeV) neutron
% c = 2200/300E6; %2200m/s
% ke = -(mass*(c^2 + (1 - c^2)^(1/2) - 1))/(c^2 - 1);