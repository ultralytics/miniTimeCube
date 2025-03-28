% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

clc; clear; close; format long g
MeV = 4;
e = 4; %event number
load(['savedEvents' num2str(MeV) 'MeV.mat'])

% E = sqrt(P^2*c^2+(m*c^2)^2)
% P = sqrt(E^2+(m*c^2)^2)/c
% E1^2 - (m1*c^2)^2  =  E2^2 - (m2*c^2)^2

mass.neutron = 939.56560; %MeV/c^2
mass.positron = .510998918; %MeV/c^2
mass.proton = 938.272013; %MeV/c2
c = 299792458; %m/s

constant    = 0.782588081999961; % mass.neutron-mass.positron-mass.proton

rpy         = [0 90 0]*pi/180;
C           = fcnRPY2DCM_B2W(rpy);

positronUnitVector  = C * GEANT.positron(e).uVec.CC(1,1:3)';
neutronUnitVector   = C * GEANT.neutron(e).uVec.CC(1,1:3)';
p1 = sqrt( (GEANT.positron(e).uVec.energy + 0/10000)^2 - mass.positron^2 ) * positronUnitVector;
p2 = sqrt( (GEANT.neutron(e).uVec.energy + 0/10000)^2 + constant^2 - mass.neutron^2  ) * neutronUnitVector; %SHOULD THIS CONSTANT BE IN HERE?

fprintf('p1=%18.10f%18.10f%18.10f\np2=%18.10f%18.10f%18.10f\np1+p2=%15.10f%18.10f%18.10f\n',p1,p2,p1+p2)

%ENERGY EQUATIONS ---------------------------------------------------------
%1) that the total energy on both sides had better be equal:
%neutrino kinetic energy (~ total energy as neutrino mass is negligible) + proton mass =
%%neutron total energy + positron total energy
Energy0 = mass.proton + MeV %VERIFIED CORRECT
Energy1 = GEANT.positron(e).uVec.energy + GEANT.neutron(e).uVec.energy %%VERIFIED CORRECT


