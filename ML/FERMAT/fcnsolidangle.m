function sa = fcnsolidangle(r,A)
%A=area of surface (units^2)
%r=radius (units)
sa = 0.5*(1 - r/sqrt(r^2+A/pi));