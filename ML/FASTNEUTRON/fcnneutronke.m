function ke = fcnneutronke(x)
%x=[xyz t E; xyz t E];
neutronmass = 939.565378; %MeV/c^2
c = 299.792458; %mm per ns, speed of light
dt = x(2,4) - x(1,4); %ns between bounce1 and bounce2
dx = fcnrange(x(2,1:3),x(1,1:3)); %mm
v = dx./dt; %mm/ns
ke = (1/2*neutronmass/c^2)*v.^2; %ke = 1/2m*v^2

