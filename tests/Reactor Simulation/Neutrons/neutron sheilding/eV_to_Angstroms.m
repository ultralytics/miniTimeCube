%Function for calculating the DeBroglie wavelength of a neutron of know Kinetic Energy 

%for use in...http://www.ncnr.nist.gov/instruments/bt1/neutron.html in
%determing the extinction of neutrons in a material

function [Debroglie_wl_Angstroms ] = eV_to_Angstroms(neutron_KE_eV)

eV_per_Joules = 6.24150974e18;%[eV/Joule]
c = 2.99792458e8;%m/s
h = 6.62606957e-34;%Joule-s
h_eV = h*eV_per_Joules;

%neutron rest mass:
E0 = 939.565; % MeV
m0 = 1.674927351e-27; %kg

%proton rest mass: (for validation in book:  good to ~.005
%E0 = 938.272; % MeV
%m0 = 1.672621777e-27;% kg

%KE = 617; % MeV

KE = neutron_KE_eV*1e-6;

%KE_per_E0 = [ 0:.2:1.4];
%KE = KE_per_E0*E0;

frac = ((KE./E0)+1).^(-2);
v_per_c = sqrt(1-frac) ;
vel = c*v_per_c;%this looks like it works. (At least it agrees with the graph in the mod. phys book.

gamma = (1-(v_per_c.^2)).^-.5;
 
wl = h./(gamma.*m0.*vel);
Debroglie_wl_Angstroms = wl*1e10;

return
