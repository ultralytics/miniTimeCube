function [] = simpleBackprojectionExample()
clc; close all; clear;
reset(RandStream.getGlobalStream)

%DEFINE CONSTANTS ---------------------------------------------------------
input.Material(1).mu(1) = 1.58; %mean index of refraction
input.scintillatorDecay = .3; %(ns) exponential scintillation delay
timeStampNoise = .1; %ns
smearExp = input.table.smearExp;
c = input.Material(1).mu(6); %group speed of light through medium (mm/ns);

%DEFINE DETECTOR PIXELS ---------------------------------------------------
np = 30; %number of pixels per face
L = 2000; %detector side length (mm)
Lh = L/2; %half length
x = linspace(-Lh+Lh/np,Lh-Lh/np,np);
y = ones(1,np);
pixels.x = [x        x        y*Lh    y*-Lh];
pixels.y = [y*Lh     y*-Lh    x       x   ];

%DEFINE VOXELS ------------------------------------------------------------
nv = 80; %number of voxels per side
x = linspace(-Lh+Lh/nv,Lh-Lh/nv,nv);
[image.x, image.y] = ndgrid(x,x); %voxel centers
image.xv = image.x(:); %x vector
image.yv = image.y(:); %y vector

%CREATE MEASUREMENTS ------------------------------------------------------
nz = 50; %number of measurements
event.t = 2; %true time event happens
event.xy = [0 0]; %true event location within detector (mm)
zid = ceil(rand(nz,1)*np*4); %measured pixel id
zt = event.t + sqrt( (event.xy(1) - pixels.x(zid)).^2 + (event.xy(2) - pixels.y(zid)).^2 )/c + randn(1,nz)*timeStampNoise; %measured time

fig;
plot(pixels.x,pixels.y,'ks');
plot(pixels.x(zid),pixels.y(zid),'ks','MarkerFaceColor','g');
plot(event.xy(1),event.xy(2),'gx','MarkerSize',15)
hs = pcolor(image.x,image.y,zeros(nv,nv)); shading flat
xlabel('x (mm)'); ylabel('y (mm)'); axis equal tight
ht = title(' '); 
legend('Pixels','Photon Hits','True Point','Backprojection'); 
colorbar; %set(gca,'clim',[0 max(smearExp.ys)*nz*.5])

for tplot = event.t % + linspace(-3,3,60)
    %BACKPROJECTION -----------------------------------------------------------
    pdf = zeros(nv^2,1);
    dt = zt - tplot; %delta-time between tplot and zt
    for i=1:nz
        j = zid(i);
        r = sqrt( (image.xv - pixels.x(j)).^2 + (image.yv - pixels.y(j)).^2 ); %range (mm) from all voxels to pixel j
        %pdf = pdf + normpdf(r,dt(i)*c,timeStampNoise*c); %normal
        pdf = pdf + interp1c(smearExp.x,smearExp.pdf, dt(i)-r/c);
    end
    pdf = reshape(pdf,[nv nv]);
    
    %PLOT ---------------------------------------------------------------------
    deleteh(hs);
    hs = pcolor(image.x,image.y,pdf); shading flat
    set(ht,'string',sprintf('t=%.3gs backprojection',tplot));
    drawnow('update')
end
end


function y = normpdf(x,mu,sigma)
    y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
end

function y = exppdf(x,lambda)
y = zeros(size(x));
i = find(x>0);
y(i) = lambda*exp(-lambda*x(i));
end