function wlpdf = cherenkov_magnitude(wl,ir,speed)
%wl (nm)
%speed (c)
%wlpf (#photons/mm^2)

%perm=0.999992; %perm - relative permeability
%q=1; %q - charge in units of electron charge
%alpha=1/137.035999074;  %unitless
%wlpdf = (1E6*(2*pi)*alpha*q^2*perm)*((wl.^-2).*(1-1./(speed*ir).^2));  %Number of Photons per path length per wavelength, Num/Length(mm)/Wavelength(nm)

wlpdf = -45850.251643*(1./(ir.^2.*speed.^2) - 1)./wl.^2;  %Num/Length(mm)/Wavelength(nm)

i = ir < 1./speed;    %Determine Wavelengths where Cherenkov is impossible
if any(i(:))
    wlpdf(i) = 0;
end


