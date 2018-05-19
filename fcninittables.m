function input = fcninittables(input)
input.cube.prettyname = regexprep(input.cube.MCPname,{'.mat','_'},{'',' '});
prettyname = regexprep(input.volumeFilename,{'.mat','_'},{'',' '});

%WAVELENGTH TABLE 
X = (82:1:2000)'; %wavelength (nm)
zv = zeros(size(X));
AL0 = 2.2; %meters (BCF10=2.2, BCF12=2.7);
RF0 = 0.85; %re-emission

%QUANTUM EFFICIENCY TABLE 
switch input.cube.MCPname
    case 'manual'
        load([input.directory '/LIBRARIES/PMTS/PLANACON_XP85012-A1.mat']); %#ok<*NODEF>
    otherwise
        load([input.directory '/LIBRARIES/PMTS/' input.cube.MCPname]);
end
qe=interp1c(x,qe,X);  %qe=qe/max(qe)*0.21;  %table max is 0.21
if isfield(input.cube,'pmt'); qe=qe*input.cube.pmt.openarea; end

np = numel(X);
switch input.volumeFilename
    case 'Manual'
        load([input.directory '/LIBRARIES/MATERIALS/EJ-254.mat']);
        otherwise
        load([input.directory '/LIBRARIES/MATERIALS/' input.volumeFilename]);
end
al = interp1c(x,al,X);  al=al/weightedMean(al,scint)*AL0; %al=al/max(al)*AL0;
re = interp1c(x,re,X);  re=re/max(re)*RF0;
scint = interp1c(x,scint,X);  scint=scint./sum(scint)/1;%normalize pdf, integration size=1nm


%ADD TO INPUT STRUCTURE
A.name = prettyname;  A.PEsensor = false;
A.yield = input.yield;
A.X=[fcnsellmeier(X,'EJ-254') al zv re scint zv]; %[1-ir, 2-al, 3-at, 4-re, 5-scintillation, 6-groupSpeed (mm/ns)]
A.mu = [];

a = buildScintillatorDecay(input.scintillatorDecay);
A.R.scintT = prepareRCDF(a.x,cumsum(a.pdf));
A.R.scintWL = prepareRCDF(X,cumsum(scint));
A.cherenkov=[];
input.Material(1)=A;

%OTHER MATERIALS
A.X=[fcnsellmeier(X,'N-BK10') al zv re scint zv];
A.name = 'Schott Glass N-BK10';  A.PEsensor = true;
input.Material(2)=A;

A.X=[fcnsellmeier(X,'Air') al zv re scint zv];
A.name = 'Air';  A.PEsensor = false;
input.Material(3)=A;

al5 = al/weightedMean(al,scint)*2.2; %(2.2 m scintillation mean)
A.X=[fcnsellmeier(X,'PS') al5 zv re scint zv];
A.name = 'PS (n=1.59)';  %http://www.crystals.saint-gobain.com/Scintillating_Fiber.aspx
input.Material(4)=A;

A.X=[fcnsellmeier(X,'PMMA') al5 zv re scint zv];
A.name = 'PMMA Cladding (n=1.49)';
input.Material(5)=A;

ir=fcnsellmeier(X,'PMMA')*1.42/1.49;
A.X=[ir al5 zv re scint zv];
A.name = 'Acrylic Outer Cladding (n=1.42)';
input.Material(6)=A;


w = scint.*qe; w = w/sum(w); %weights
for i=1:numel(input.Material)
    A = input.Material(i);
    
    ir = A.X(:,1);
    gs = groupSpeed(X,ir);
    al = A.X(:,2); %(m)
    at = al./gs*1000; %(ns), mean attenuation time = attenuationLength in meters*1000 / photon velocity in mm/s
    A.X(:,3) = at;
    A.X(:,6) = gs;
    
    A.mu = w'*A.X;
    A.cherenkov = cherenkovTable(X,ir);
    
    input.Material(i)=A;
end

input.table.smearExp = a;

input.wl = X;
input.colormap = interp1(1:np, fcnWL2RGB(X), linspace(1,np,1024));
input.cube.QE = qe;
input.cube.QEmean = sum(scint.*qe/sum(scint));
end

function speed=groupSpeed(wl,ir)
%https://en.wikipedia.org/wiki/Group_velocity
%wl = wavelength (nm)
c = 299.792458; % (mm/ns)
speed=c./(ir-wl.*[0; diff(ir)]);  %Photon Group Speed (mm/ns)
end

function c=cherenkovTable(X,ir) %CHERENKOV PHOTON DENSITY TABLE (#/mm per speed and wavelength) -----------
c.speed = linspace(0,1,101); %(c)
c.wlpdf = cherenkov_magnitude(X, ir, c.speed); %number of photons/mm
c.wlpdf(1:81,:)=0;
c.wlcdf = cumsum(c.wlpdf);
c.npe = c.wlcdf(end,:);
end

function A = buildScintillatorDecay(tf)
np = 10000;
%tf = 2.2; %2.2ns fall time
tr = 0.9;  % 0.9ns rise time
s  = 0.100;  % 1sigma smear (ns)

if tf==0 %special test case, zero delay
    x = linspace(-1,1,np)';  [~,i]=min(abs(x));
    y = zeros(size(x));  y(i)=1;
    ys = fcnpdfnorm(x,x',s)*y;  %smeared
else
    format = 'DoubleExponential'; %'RiseFall'
    switch format
        case 'EMG'
            x = linspace(-8*s-1,12*tf+3*s+1,np)'; 
            y = fcnpdfEMG(x,0,0,1/tf);
            ys = fcnpdfEMG(x,0,s,1/tf); %smeared
        case 'DoubleExponential'
            x = linspace(-10*s-25,max(12*tf,tr*6)+50,np)';
            %y = exp(-x/tf)-exp(-x/tr);  %old way
            y = exp(-x/tf).*(1-exp(-x/tr)); %GEANT way
            y(x<0)=0;
            ys = fcnpdfnorm(x,x',s)*y;  %smeared
    end
end
dx=x(2)-x(1);
y=y/sum(y)/dx;
ys=ys/sum(ys)/dx;

if all(y==0 | isnan(y)); [~,i]=min(abs(x)); y=zeros(np,1); y(i+[-1 0 1])=1;  y=y/sum(y)/dx; ys=y; end

i = ys>0;  ys(~i)=min(ys(i));
A.x = x;
A.pdf = y;
A.pdfsmeared = ys;

%CHERENKOV TIME 
x = linspace(-5,5,3000); %ns
y = fcnpdfnorm(x,0,.02);  i=y>0;
A.cx = x(i);
A.cpdf = y(i);

%CHERENKOV ANGLE
x = linspace(-1,1,3000); %ct
y = fcnpdfnorm(x,.6,.2);  i=y>0;
A.ax = x(i);
A.apdf = y(i);
end

