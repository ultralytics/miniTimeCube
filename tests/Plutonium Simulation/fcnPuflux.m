function [Bi, Si, Ci, Cisb] = fcnPuflux()
t =     60; % (s)
n =     61000; % (n/s/kgWGP) Luminosity Weapons Grade Plutonium WGP = 6% Pu240, 94% Pu239
area =  13^2; %detector area (cm^2)
es =    .16; %signal efficiency
eb =    .065; %background efficiency
mass =  .0005; %kg WGP
range = 30; %cm
al = 215*100; %(cm) attenuation length

Si = mass*n*exp(-range/al)./(4*pi*( range ).^2)*area*t*es; %signal neutrons (n/area/time)
Bi = .0134*(area*3)*t*eb; %.0134 Goldhagen (n/cm^2/s) %Only multiply by 3 sides per Pieter Mumm
Bi = Bi/3; %SATO OVER WATER https://mail.google.com/mail/u/0/?shva=1#search/water/14056204ed8c753b

Cisb(1) = fcnRandomPoisson(Si);
Cisb(2) = fcnRandomPoisson(Bi);
Ci = sum(Cisb);
end

