function k = fcnoptimizerk(input, output)
%MTCflag = ischecked(evalin('base','handles.GUI.realdataflag'));
MTC766flag = false;
if MTC766flag
    k.QEmap = input.cube.QEmean * ones(input.cube.pixels,1);  [~,i]=fcnpruninglist;  k.QEmap(~i) = 0;
else
    k.QEmap = input.cube.QEmean * ones(input.cube.pixels,1);
end
j = find(k.QEmap>0);  k.QEmap=k.QEmap(j);  k.pid=j;
k.npixels = numel(j);

k.fi = input.cube.all.fi(j); %face index
k.nxyz = -input.cube.all.normalVec(j,:);
k.pxyz = input.cube.all.xyz(j,:);

k.fn = cell(6,1);  k.fp=k.fn;
for i=1:6
    k.fn{i} = k.nxyz(k.fi==i,:); %face normals
    k.fp{i} = k.pxyz(k.fi==i,:); %face pixel pos
end

k.smearExp = input.table.smearExp;
k.al = input.Material(1).mu(2)*1000; %(mm now from meters)
k.nial = -1/k.al; %negative inverse attenuation length (attenuation length in m)
k.pa = input.cube.pixelarea; %pixel area
k.paofp = k.pa/(4*pi); %pixel area over 4pi
k.res  = k.pa/pi; %pixel area equivalent radius r squared

k.qemean = input.cube.QEmean;
k.yield = input.Material(1).yield;
k.ci  = 1/input.Material(1).mu(6); %1/groupSpeed 

if nargin==2
    t=output.ZI(j,4);  t(t==0)=nan;
    [~,inliers] = fcnsigmarejection(t); t(~inliers)=0;  k.zt=t;
    k.ztpid = find(k.zt~=0);
    
    k.nhpp = output.N(j); %number of hits per unique pixel
    k.poiss.gammaln = gammaln(k.nhpp+1);
else
    k.nhpp = ones(1,input.cube.pixels);
end

k.ior1 = input.Material(1).mu(1);
k.ior2 = input.Material(2).mu(1);

k.ctx = linspace(0,1,3000);
k.yFr = fcnfresnelc(k.ior1,k.ior2,k.ctx)';
k.yFt = 1 - k.yFr;

k.shapeID = input.cube.shapeID;
k.Lr = input.cube.Lr;
if input.cube.shapeID==3 %cylinder
    k.radius = input.cube.radius;
    k.height = input.cube.height;
end

k.timeflag = false;
k.reflections = 0; %number of reflections to model
k.cherenkov = 0;