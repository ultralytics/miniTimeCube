function e = fcnwl2e(wl,ir) %wavelength (nm), index of refraction
c = 299792458; %m/s, speed of light
planck = 4.13566751691E-15; %eV*s
e = (planck*c*1E3)./(wl.*ir); %MeV, includes nm2m conversion
end

