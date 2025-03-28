% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [C, XC, TC, EC, EFC, str] = fcngetcandidates(input,N,X,T,F,particleName)
%C=candidates, XC = candidate estimates, TC = candidate truths, EC = candidate estimate errors (xhat-xtrue)
%N=number of photons, X=xhat estimate, T=xtrue, E=error (xhat-xtrue), F=estimator failure flag, str='antineutrino'
nmc = size(X,1); %number of mc runs
nxp = size(X,2); %number of xhat parameters
nts = size(X,3); %number of ts runs

switch particleName
    case {'antineutrino','Antineutrino'}
        wc = 5; %(mm) wall cut [min]
        pndc = [0 1000]; %(mm) positron-neutron distance [min max]
        pec = [1.000 8.000]; %(MeV) prompt energy [min max]
        dec = [0.040 0.400]; %(MeV) delayed energy [min max]
        dtc = [50 12000]; %(ns) delay time [min max]
        ppc = [20 1E4]; %prompt photon count
        dpc = [20 400]; %delayed photon count
        
        %FOR LITHIUM DOPING
        dec = [0.200 0.600]; %(MeV) delayed energy [min max]
        dpc = [20 900]; %delayed photon count
        
        ncc = 8; %number of candidate criteria
        C = false(nmc,ncc,nts);  XC.mu = zeros(nts,nxp);  XC.sigma=XC.mu;  TC=XC;  EC=XC;  EFC=XC;
        for i = 1:nts
            Ep = X(:,4+0,i); %(MeV) energy of prompt
            Ed = X(:,4+5,i); %(MeV) energy of delayed
            xp = X(:,(1:3)+0,i); %(mm) position of prompt
            xd = X(:,(1:3)+5,i); %(mm) position of delayed
            dt = X(:,10,i)-X(:,5,i); %(ns) delayed-prompt time gap
            
            C(:,1,i) = Ep>pec(1) & Ep<pec(2); %prompt energy
            C(:,2,i) = Ed>dec(1) & Ed<dec(2); %delayed energy
            C(:,3,i) = fcninsidedetector(xp, input, wc) & fcninsidedetector(xd, input, wc); %wall distance cuts
            C(:,4,i) = dt>dtc(1) & dt<dtc(2); %(ns) prompt-delayed delta time
            C(:,5,i) = rangec(xd-xp)>pndc(1) & rangec(xd-xp)<pndc(2); %delayed-prompt distance
            C(:,6,i) = N(:,1,i)>ppc(1) & N(:,1,i)<ppc(2); %number of prompt photons
            C(:,7,i) = N(:,2,i)>dpc(1) & N(:,2,i)<dpc(2); %number of delayed photons
            
            [XC,TC,EC,EFC,C] = cstats(XC,TC,EC,EFC,C,X,T,i);
        end
        
        str = {sprintf('E_p:       %.3g-%.3gMeV',pec),...
            sprintf('E_d:       %.3g-%.3gMeV',dec),...
            sprintf('wall:      %.3gmm',wc),...
            sprintf('dt:        %.3g-%.5gus',dtc/1E3),...
            sprintf('dx:        %.3g-%.3gmm',pndc),...
            sprintf('PE_p:      %.3g-%.5g',ppc),...
            sprintf('PE_d:      %.3g-%.5g',dpc),...
            sprintf('combined')};
    case {'FiberGamma','Gamma','gamma'}
        
    case {'FiberNeutron','Neutron','neutron'}
        ncc = 13; %number of candidate criteria
        C = false(nmc,ncc,nts);  XC.mu = zeros(nts,nxp);  XC.sigma=XC.mu;  TC=XC;  EC=XC;  EFC=XC;
        
        ec = 2; %edge cut (mm)
        for i = 1:nts
            r=fcnrange(X(:,1:3,i)-X(:,6:8,i));  dt=X(:,9,i)-X(:,4,i);  vel=r./dt;
            C(:,1,i) = ~all(X==0,2) & ~any(isnan(X),2); %single scatter
            C(:,2,i) = X(:,14,i)>1; %> number recoils
            C(:,3,i) = all(abs(X(:,1:3,i))<(input.cube.Lr(3)-ec),2) & all(abs(X(:,6:8,i))<(input.cube.Lr(3)-ec),2); %wall cut (mm)
            C(:,4,i) = X(:,11,i)./X(:,12,i)>.1; % PE ratio
            C(:,5,i) = X(:,5,i)>.01 & X(:,5,i)<10; % E1 (MeV)
            C(:,6,i) = X(:,10,i)>.01 & X(:,10,i)<10; % E2 (MeV)
            C(:,7,i) = X(:,11,i) > 10; % P1 PE count (was 50)
            C(:,8,i) = X(:,12,i) > 10; % P2 PE count (was 50)
            C(:,9,i) = r > 30; %dx (mm)
            C(:,10,i) = dt>.1 & dt<20; %dt (ns) (was 1-20)
            C(:,11,i) = vel<40 & vel>5; %V1 (mm/ns)
            C(:,12,i) = X(:,17,i)<6 & X(:,17,i)>0; %TRUE E
            str='';
            [XC,TC,EC,EFC,C] = cstats(XC,TC,EC,EFC,C,X,T,F,i);
        end
    case {'MTCneutron','MTCNeutron'}
        ncc = 13; %number of candidate criteria
        C = false(nmc,ncc,nts);  XC.mu = zeros(nts,nxp);  XC.sigma=XC.mu;  TC=XC;  EC=XC;  EFC=XC;
        
        ec = 2; %edge cut (mm)
        for i = 1:nts
            r=fcnrange(X(:,1:3,i)-X(:,6:8,i));  dt=X(:,9,i)-X(:,4,i);  vel=r./dt;
            C(:,1,i) = ~all(X==0,2) & ~any(isnan(X),2); %single scatter
            C(:,2,i) = X(:,14,i)>1 & X(:,14,i)<6; %> number recoils
            C(:,3,i) = all(abs(X(:,1:3,i))<(input.cube.Lr(3)-ec),2) & all(abs(X(:,6:8,i))<(input.cube.Lr(3)-ec),2); %wall cut (mm)
            C(:,4,i) = X(:,11,i)./X(:,12,i)>.1; % PE ratio
            C(:,5,i) = X(:,5,i)>.2 & X(:,5,i)<4; % E1 (MeV)
            C(:,6,i) = X(:,10,i)>.2 & X(:,10,i)<4; % E2 (MeV)
            C(:,7,i) = X(:,11,i) > 25; % P1 PE count (was 50)
            C(:,8,i) = X(:,12,i) > 15; % P2 PE count (was 50)
            C(:,9,i) = r > 10; %dx (mm)
            C(:,10,i) = dt>1 & dt<20; %dt (ns) (was 1-20)
            C(:,11,i) = vel<35 & vel>2; %V1 (mm/ns)
            C(:,12,i) = X(:,17,i)<10 & X(:,17,i)>.5; %E0hat
            str='';
            [XC,TC,EC,EFC,C] = cstats(XC,TC,EC,EFC,C,X,T,F,i);
        end
end
end

function [XC,TC,EC,EFC,C] = cstats(XC,TC,EC,EFC,C,X,T,F,i)
j = all(C(:,1:end-1,i),2); if ~isempty(F); j=j & ~F(:,1,i) & ~F(:,2,i); end %all cuts & prompt & delayed convergence;
C(:,end,i) = j; %overall candidates

Xji = X(j,:,i);
Tji = T(j,:,i);

%XHAT CANDIDATE STATS
XC.mu(i,:) = mean(Xji,1);
XC.sigma(i,:) = std(Xji,[],1);

%XTRUE CANDIDATE STATS
TC.mu(i,:) = mean(Tji,1);
TC.sigma(i,:) = std(Tji,[],1);

%ERROR (XHAT-XTRUE) CANDIDATE STATS
e = (Xji - Tji); %error
EC.mu(i,:) = mean(e,1);
EC.sigma(i,:) = std(e,[],1);

%ERROR FRACTION (XHAT-XTRUE)/XTRUE CANDIDATE STATS
bias = EC.mu(i,:);
e = abs(((Xji - Tji)-bias)./Tji); %error
EFC.mu(i,:) = mean(e,1);
EFC.sigma(i,:) = std(e,[],1);
end