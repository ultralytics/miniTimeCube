function wl = fcne2wl(e,ir) %energy (MeV)
c = 299792458; %m/s, speed of light
planck = 4.13566751691E-15; %eV*s
wl = (planck*c*1E3)./(e.*ir); %nm, includes nm2m conversion
end

