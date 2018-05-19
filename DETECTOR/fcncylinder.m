function [P, input] = fcncylinder(input)
P = cell(1,6);
d2r = pi/180;
input.cube.requestedpixels = max(input.cube.requestedpixels,100);
r = input.cube.Lr(1); input.cube.radius = r;
h = input.cube.Lr(3)*2; input.cube.height = h;

input.cube.volume = pi*(r/1000)^2*h/1000;
sidefacewidth = 2*r*pi/4;
input.cube.area.perface = [ sidefacewidth*h*[1 1 1 1]  pi*r^2*[1 1] ];
input.cube.area.all = sum(input.cube.area.perface);

switch input.cube.MCPname
    case 'manual'
        input.cube.pixelsize = sqrt(input.cube.area.all*input.cube.coverageFraction/input.cube.requestedpixels) * [1 1];
        input.cube.footprint = input.cube.pixelsize;
    otherwise
end
ppf = round(input.cube.area.perface / prod(input.cube.footprint) * input.cube.coverageFraction); %pixels per face

w2h = sidefacewidth/h;
ppc = round(sqrt(ppf(1)/w2h)); %pix per col
ppr = round(ppf(1)/ppc); %pix per row

azdegpp = atan(input.cube.footprint(1)/2/r)*r2d; %az deg per tile
ppr = min(ppr,floor(45/azdegpp));
ppc = min(ppc,floor(h/(input.cube.footprint(2))));

az0vec = [0 90 180 -90];   
for i = 1:4
    az = (midspace(-45,45,ppr)'+az0vec(i))*d2r;
    z = midspace(-1,1, ppc)'*h/2;
    
    el = atan(z/r);
    
    [EL, AZ] = ndgrid(el, az);
    R = repmat(sqrt(r^2+z.^2),[1 ppr]);
    P{i} = fcnSC2CC(R(:),EL(:),AZ(:));
end

%TOP AND BOTTOM
ppr = round(sqrt(ppf(5)/0.7854))-1;
x = midspace(-1,1,ppr)*(r-input.cube.footprint(1)/4);
[X, Y] = ndgrid(x,x);
R = sqrt(X(:).^2+Y(:).^2);
i = find(R < (r-input.cube.footprint(1)/sqrt(2)));
np = numel(i);

P{5} = [X(i) Y(i) repmat(-h/2,np,1)]; %5-TOP
P{6} = [X(i) Y(i) repmat( h/2,np,1)]; %6-BOTTOM


input.cube.area.perface = input.cube.area.perface/1000^2; %mm^2 to m^2
input.cube.area.all = input.cube.area.all/1000^2;
end

