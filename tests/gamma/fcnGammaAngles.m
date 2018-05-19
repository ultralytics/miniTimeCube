function [] = fcnGammaAngles(input,flags,handles,GEANT,G1)
clc
%fig; plot(theta*r2d,de,'.','markersize',25)

%A=zeros(10000,6);
A=cell(50000,1);
parfor j=1:size(A,1);
    j
    [input2,output2,G1,photons2,PE,flags2] = fcn1Event(input,flags,handles,GEANT,[0 0 1 j]);
    
    i=G1.tid==1 & G1.pid==22; %gammas
    nv = find(G1.tid==1);  if numel(nv)<3; fprintf('insufficient GEANT TIDs.\n'); continue; end %neutron rows

    de=[G1.ke1(i)-G1.ke2(i)]; %MeV
    v=[G1.p2(i,:)-G1.p1(i,:)];

    theta=fcnangle(v(1:end-1,:),v(2:end,:));  de=de(1:end-1);  ke1=G1.ke1(i);  ke1=ke1(1:end-1);
    
    A{j}=[ke1 de theta];
    
    
%     i=G1.p2inside & G1.parentinside & G1.pid==22;
%     if numel(de)==0; theta=0; de=0; end
%     if PE.n==0; npixels=0; else npixels=numel(fcnunique(PE.pixel)); end
%     B=[de(1) theta(1) double(PE.n), npixels, sum(G1.ke1(i)-G1.ke2(i)), sum(double(i))];
%     A(j,:) = B;
end
%B=cat(1,A{:}); save gamma662keVStats.mat B
A=[]; 
load('gamma3MeVStats.mat'); A=[A; B];
load('gamma10MeVStats.mat'); A=[A; B];
load('gamma662keVStats.mat'); A=[A; B]; B=A; clear A
i=B(:,1)==.662; %B(:,1)>.4 & B(:,1)<.405;
fig; histogram(B(i,2),100)

[N,xe,ye]=histcounts2(B(i,1),B(i,2),200); x=xe(2:end)/2+xe(1:end-1)/2; y=ye(2:end)/2+ye(1:end-1)/2;
fig; imagesc(x,y,N'); fcntight('c'); xyzlabel('E_\gamma (MeV)','E_{e^-} (MeV)')

% H=fspecial('gaussian',3,1); 
% N=imfilter(N,H,'replicate');
% fig; imagesc(x,y,N'); fcntight('c'); xyzlabel('E_\gamma (MeV)','E_{e^-} (MeV)')

F=griddedInterpolant({x,y},N,'linear','none');
xi=linspace(0,6,1000)'; 
zi=F(xi,zeros(size(xi))+.3); fig; plot(xi,zi); xyzlabel('E_\gamma (MeV)')

fig; 
for a=[.5 1 1.5 2]
    plot(xi,F(zeros(size(xi))+a,xi),'DisplayName',sprintf('%.3f MeV \\gamma',a));
end
xyzlabel('E_{e^-} (MeV)'); legend show; fcnlinewidth(2); fcnfontsize(12);


save gammaF.mat F

% fig(1,3,1);
% sca; plot(A(:,2),A(:,1),'.','markersize',25); xyzlabel('gamma d\theta','gamma dE (MeV)'); set(gca,'Xlim',[0 pi]);
% sca; fcnhistc(cos(A(:,2))); xyzlabel('cos(\theta)'); legend('0.2 MeV \gamma first scatter angle')
% sca; fcnhistc(A(:,1)); xyzlabel('dE (MeV)'); legend('0.2 MeV \gamma first scatter dE')

% fig(2,2,'19cm'); i=A(:,6)>0;
% sca; fcnhistc(A(i,3)); xyzlabel('GEANT Co60-MTC PEs')
% sca; fcnhistc(A(i,4)); xyzlabel('GEANT Co60-MTC lit pixels')
% sca; fcnhistc(A(i,5)); xyzlabel('GEANT Co60-MTC dE (MeV)')
% sca; fcnhistc(A(:,6)); xyzlabel('GEANT Co60-MTC recoils'); mean(A(:,6)==0)


% % %COMPTON SCATTER EQUATION -------------------------------------------------
% % %https://en.wikipedia.org/wiki/Compton_scattering
% c = 299792458; %m/s, speed of light
% h = 4.13566751691E-15; %eV*s %planck constant
% me = 0.510998910; %electron rest mass (MeV) or (MeV/c^2)
% 
% E1 = 0.500; %Energy BEFORE scatter
% E2 = linspace(.2,.500,100); %energy AFTER scatter
% 
% syms h c E2 E1 me theta
% %solve(h*c/E2 - h*c/E1 - h/me/c*(1-cos(theta)),theta)
% 
% ct = me*(1./E1 - 1./E2) + 1;
% plot(acos(ct),E2)