function []=fcnReactorCollection()
clc; clear; close all; tic
t = 60*60*24; %s
L = 13; %Length of cube side (cm)
rp = 20; %Reactor Power (MWth)
range = 4; %Reactor-Detector range (m)
twindow = 12E-6; %collect window (seconds)
nprotons = 5.6573e+22 * L^3; %5.6573e+22 protons per cm^3
shield.lead = 0; %cm
shield.polyethylene = 0; %cm
shield.MWE = 0; %m
fprintf('Starting Reactor-Detector Simulation ...\nAssuming %.3gcm detector, %.3gL volume containing %.2e protons\nAssuming %.3gcm lead, %.3gcm polyethelyne and %.3gMWE shielding\nCollecting for %.4gs at %.3gm from a %.3gMWth reactor\n\n', ...
    L,L^3/1E3,nprotons,shield.lead,shield.polyethylene,shield.MWE,t,range,rp)

s2days = 1/60/60/24;
area = L^2;
np = zeros(6,7);

s = fcngetobservedspectra(shield,true);
ns = numel(s);
flux = zeros(ns,1);

snuebar = fcnspec1s(1.8:.01:11, range/1000, [], rp/1000, nprotons, 1*s2days);

%fcnspec1s(1.8:.01:11, .004, [], .020, 1.243E26, 1)
%ia = 63.5 * 66.04 * 119.38 / 1E3;  %%Inner Air Volume of Cave in Liters
%iar = 2.2/500; %0.0044

flux(1) = snuebar.n; %Antineutrinos/s
flux(2) = 1.3 * area * rp/20; %Reactor Neutrons/s
flux(3) = 0.7 * area; %Atmospheric Neutrons/s
flux(4) = 51 * area * rp/20; %Reactor Gammas/s
flux(5) = 25 * area; %Atmospheric Gammas/s
flux(6) = .0165 * area; %Muons/s  .0165/cm^2/s
%flux = flux.*(1-[0 .9869 .9869 .9923 .9923 .6242]'); %Shielding fractions


%flux = [snuebar.n, 0, 0.05, 0, 13, 4]';
flux = [snuebar.n/area/6, 0, 9.7E-5, 0, 1.1E-2, 2E-3] * area*6;


flux = flux*t;
fprintf('Detector Flux in %.4gs:\n',t)
for i=1:6
   fprintf('%s: %.4g\n',s(i).simplename,flux(i)) 
end

zv = zeros(1E7,1,'uint8');
particleID = zv;
candidateID = zv;
%particleE = zv;
clear zv

ne = 0;
for i=1:6 %particle
    cf = s(i).sf/100;
    fluxi = flux(i)*cf;
    %pdfx = s(i).pdfx;  dx = pdfx(2)-pdfx(1);
    for j=1:7 %candidate type 1. Full, 2. Prompt->Prompt, 3. Prompt->Delayed, 4. Delayed->Prompt, 5. Delayed->Delayed, 6. Unobserved, 7. Rejected, 8. Shielded
        np(i,j) = fcnRandomPoisson(fluxi(j));
        if np(i,j)>0
            %cdfy = cumsum(s(i).pdfyc(:,j))*dx;
            %particleE(k,1) = fcnrandcdf(cdfy, pdfx, np(i,j));
            k = ne+1:ne+np(i,j);
            particleID(k,1) = i;
            candidateID(k,1) = j;
            ne=ne+np(i,j);
        end
    end
end
T = rand(ne,1)*t;
particleID = particleID(1:ne);
candidateID = candidateID(1:ne);

[~, i] = sort(T); T=T(i); particleID=particleID(i); candidateID=candidateID(i); clear i
fprintf('\n%.5g total events in detector in %.10gs ...\n',ne,t)

i = candidateID~=6;  n=sum(i);  T=T(i);  particleID=particleID(i); candidateID=candidateID(i); %observed events
fprintf('%.5g observed events (<%.2g%% duty cycle) ...\n',n,n*twindow/t*100)

dtiw = diff(T)<twindow; %dt in window

l=0;
y=zeros(1E7,1);
flags = true(n,1); nkv=zeros(n,1,'uint8');
for i=1:n-1
    if flags(i) && dtiw(i)
            j = (i+1):(i+min(100,n-i));
            k = (T(j)-T(i))<twindow; nk=sum(k);
            
%                 if flags(i) && (T(i+1)-T(i)) < twindow
%             j = (i+1):(i+min(100,n-i));
%             k = (T(j)-T(i))<twindow; nk=sum(k);
            
            nkv(i)=nk;
            
            if nk==1
                l=l+1;
                y(l) = i;
                flags(i+1) = false; %skip
            else % nk>1
                flags(j(k)) = false;
            end
    end
end
y=y(1:l); nkv=nkv(nkv~=0);  n=l;
fprintf('%.4g double-bookends in %.6gus windows ...\n',n,twindow*1E6)

i = candidateID(y)~=7 & candidateID(y+1)~=7;
fprintf('%.4g candidate (unrejected) double-bookends ...\n',sum(i)); y=y(i);

cid = candidateID(y);  i = find(cid==2 | cid==4);  n=numel(i);  y=y(i);
fprintf('%.4g double-bookends have prompt candidates first ...\n',n)

cid = candidateID(y+1);  i = find(cid==3 | cid==5);  n=numel(i);  y=y(i);
fprintf('%.4g double-bookends have delayed candidates second ...\n\nThese final candidates are composed of the following particle pairs:\n',n)


a = particleID(y); %prompt particles
b = particleID(y+1); %prompt particles

c = unique([a b],'rows');
for i=1:size(c,1)
   nc = numel(find(a==c(i,1) & b==c(i,2)));
   fprintf('%7g %25.25s - %.25s\n',nc, s(c(i,1)).simplename(1:end-1), s(c(i,2)).simplename(1:end-1))
end

i=find(candidateID==1);  m=numel(i);
fprintf('\nIn addition the %.6g uncorrelated events above, the following %.6g correlated events were observed:\n',n,m)
i=find(particleID==1 & candidateID==1);  n=numel(i);
fprintf('%.4g Full-Event Antineutrino Events Observed ...\n',n)
i=find(particleID==2 & candidateID==1);  n=numel(i);
fprintf('%.4g Full-Event Reactor Neutron Events Observed ...\n',n)
i=find(particleID==3 & candidateID==1);  n=numel(i);
fprintf('%.4g Full-Event Atmospheric Neutron Events Observed ...\n',n)

fprintf('\n')
toc
end