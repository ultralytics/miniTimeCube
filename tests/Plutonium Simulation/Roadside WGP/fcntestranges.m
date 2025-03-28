% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcntestranges()
clc

nd = 30;
vd = linspace(3,30,nd);

n = 1E6;
xr = [linspace(0,0,n)' linspace(-500,500,n)' linspace(0,0,n)'];

fig; fcnplot3(xr,'r.'); zv = zeros(1,nd); fcnplot3([vd' zv' zv'],'b.'); axis equal vis3d; box on; grid on; legend('car positions','detector positions'); fcn3label('x (m)','y (m)','z (m)')

for i = 1:nd
    xd = [-vd(i) 0 0];
    dx = xr-xd;
    dxu = fcnvec2uvec(dx);
    r = fcnrange(dx);
    
    normald = repmat([1 0 0],[n 1]);
    [~, ct] = fcnangle(normald,dxu);
    sflux(i) = sum(abs(ct)./r.^2); %sum flux
end
fig; plot(vd,sflux,'r');


for i = 1:nd
    xd = [-vd(i) 0 0];
    dx = xr-xd;
    dxu = fcnvec2uvec(dx);
    r = fcnrange(dx);
    
    normald = repmat([0 1 0],[n 1]);
    [~, ct] = fcnangle(normald,dxu);
    sflux(i) = sum(abs(ct)./r.^2); %sum flux
end
plot(vd,sflux,'g');


for i = 1:nd
    xd = [-vd(i) 0 0];
    dx = xr-xd;
    dxu = fcnvec2uvec(dx);
    r = fcnrange(dx);
    sflux(i) = sum(1./r.^2); %sum flux
end
plot(vd,sflux,'b');

legend('1m^2 facing toward the road','1m^2 facing down the road','1m^2 cylindrical')
fcn3label('distance from road (m)','integrated flux')
grid on
end

