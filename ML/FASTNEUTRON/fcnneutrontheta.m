% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function theta = fcnneutrontheta(ke0,e)
ep = e(1); %energy of first proton hit
en = ke0; %energy of neutron before first proton hit
theta = asin(sqrt(ep/en))*180/pi;
