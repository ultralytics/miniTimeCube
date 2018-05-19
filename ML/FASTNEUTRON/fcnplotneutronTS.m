function [] = fcnplotneutronTS(input,MC,tsv)
if nargin==1
      [fname, pname] = uigetfile('*.mat','Pick a TS file'); 
      if isempty(fname); return; else load([pname fname]); end
end

close all; clc
[h,hf1]=fig(5,2);
c = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 1 0],[1 0 1]};
nb = 50;
x1 = linspace(-1+1/nb,1-1/nb,nb);
av = .5;
ni = numel(tsv);

name = cell(1,ni);
for i=1:ni
    name{i} = sprintf('%.0fMeV',tsv(i));
end
name{1} = 'Bk'; %FOR SPECIAL BACKGROUNDS ONLY!  
name{2} = 'Pu'; %FOR SPECIAL BACKGROUNDS ONLY!

str = cell(10,ni);
for i=1:ni
    if i>6
        c{i} = rand(1,3);
    end
    ce = .7*(1-c{i})+c{i}; %edgecolor
    a = zeros(size(MC.xhat(:,1,1),1),13);
    a(:,1) = MC.xhat(:,1,i)~=0  &  MC.xtrue(:,18,i)>0;
    a(:,2) = MC.xhat(:,6,i)>10 & MC.xhat(:,6,i)<300; %P1 > 10 photons
    a(:,3) = MC.xhat(:,12,i)>5; %P2 > 2  photons
    a(:,4) = MC.xhat(:,10,i) - MC.xhat(:,4,i) > 1; %dt > 1ns
    a(:,5) = fcnrange(MC.xhat(:,1:3,i) - MC.xhat(:,7:9,i)) > 10; %dx > 10mm
    a(:,6) = MC.xhat(:,5,i)>0.200; %.3 %first bounce > 100keV
    a(:,7) = MC.xhat(:,6,i)./MC.xhat(:,12,i) > .2; %dE1/dE2 must be greather than 0.20
    a(:,8) = MC.xhat(:,5,i)./MC.xhat(:,18,i) < .9; %dE1./E0 < 0.9
    a(:,9) = MC.xhat(:,5,i)./MC.xhat(:,18,i) > .1; %dE1./E0 > 0.1
    %a(:,10) = fcninsidedetector(MC.xhat(:,1:3,i), input, 10); %P1 wall cut
    %a(:,11) = fcninsidedetector(MC.xhat(:,7:9,i), input, 10); %P1 wall cut
    a(:,12) = MC.collectedPhotonCount(:,2,i)>2; %delayed bookend
    %a(:,13) = all(MC.protonprotonflag(:,:,i),2); %proton-proton bounces only!
    j = all(a(:,[1:9]),2);    mean(j)
    nj = sum(j);
    
    mean(a)

    
    
    MC.xhat(j,19,i) = input.neutron.dEfit(MC.xhat(j,19,i));
    p1hat = MC.xhat(j,1:3,i);
    p2hat = MC.xhat(j,7:9,i);
    pnhat = MC.xhat(j,15:17,i);
    anglehat = MC.xhat(j,13,i);
    E0hat = MC.xhat(j,18,i);
    
    p1 = MC.xtrue(j,1:3,i);
    p2 = MC.xtrue(j,7:9,i);
    pn = MC.xtrue(j,15:17,i);
    angle = MC.xtrue(j,13,i);
    E0 = MC.xtrue(j,18,i);
    
    a=MC.xhat(j,:,i) - MC.xtrue(j,:,i);
    %results.xhat  = [p1 t1 e1 nphotons1, 7p2 t2 e2 nphotons2, 13anglehat, nphat, neutronCapturexyzhat, 18E0hat, E0hatMLpoint, E1hat ];
    %results.xtrue = [p1 t1 e1 nphotons1, 7p2 t2 e2 nphotons2, 13angle,    np   , neutronCapturexyz     18E0     E0          , E1    ];
    
    if nj==0
        for si=1:size(str,1); str{si,i}=''; end
        fprintf('\nWARNING: NO CANDIDATES FOUND. NOT PLOTTING RESULTS.\n');
        continue
    end
    
    
    y(1,i) = sum(j)/numel(j); %efficiency
    
    b = MC.minangleerror(j,1,i);
    y(2,i) = fcnstdnonlin(b,[],0); %angle cone resolution
    
    b=a(:,5)./mean(MC.xtrue(j,5,i)); d=fcnsigmarejection(b,6,3);
    y(3,i) = fcnstdnonlin(b,[],0); %Neutron First Bounce Energy Resolution

    b=a(:,4); d=fcnsigmarejection(b,3,3);
    y(4,i) = fcnstdnonlin(b,[],0); %Neutron First Bounce Time Resolution
    
    b=a(:,1:3); d=fcnsigmarejection(b(:),3,3);
    y(5,i) = fcnstdnonlin(b,[],0); %Neutron First Bounce Position Resolution

    b=a(:,7:9); d=fcnsigmarejection(b(:),3,3);
    y(6,i) = fcnstdnonlin(b,[],0); %Neutron Second Bounce Position Resolution

    b=a(:,10); d=fcnsigmarejection(b,3,3);
    y(7,i) = fcnstdnonlin(b,[],0); %Neutron Second Bounce Time Resolution

    %sca(h(9))
    %b=a(:,18)./MC.xtrue(j,18,i); d=fcnsigmarejection(b,2,3); if ~exist('x7','var'); x7=linspace(min3(d),max3(d),nb); end
    %y = hist(d,x7); bar(x7,y,1,'facecolor',c{i},'edgecolor',ce); alpha(av)
    %str{9,i} = sprintf('%s: %.2f 1\\sigma',name{i},fcnstdnonlin(b,[],median(b)));
end
y = y';

b = y(:,1).*y(:,2:end).^-2;  b(~isfinite(b))=0;  %info = efficiency / inverse of variance
c = b./max(b,[],1); %normalize
str = {'angle cone','E_1','T_1','P_1','P_2','T_2'};
fig; plot(tsv,c,'.-'); legend(str); set(gca,'xtick',tsv);  xyzlabel('Trade Study Scenario','Relative Information','',MC.FileName)


% sca(h(1))
% xlabel('cos(\Theta)'); ylabel('number'); title(sprintf('Neutron Angle Errors'));
% legend(str(1,1:ni));
% sca(h(2))
% title('Neutron First Bounce {\bfEnergy Resolution}'); xlabel('energy error (fraction)'); ylabel('number'); hl=legend(str(2,1:ni)); %fcnlegend(hl,av)
% sca(h(3))
% title('Neutron First Bounce {\bfTime Resolution}'); xlabel('time error (ns)'); ylabel('number'); legend(str(3,1:ni));
% sca(h(4))
% title('Neutron First Bounce {\bfPosition Resolution}'); xlabel('position error (mm)'); ylabel('number'); legend(str(4,1:ni));
% sca(h(5))
% title('Neutron Second Bounce {\bfTime Resolution}'); xlabel('time error (ns)'); ylabel('number'); legend(str(5,1:ni));
% sca(h(6))
% title('Neutron Second Bounce {\bfPosition Resolution}'); xlabel('position error (mm)'); ylabel('number'); legend(str(6,1:ni));
% sca(h(9))
% title('E0 = dE_0^1 + E1 {\bfEnergy Resolution}'); xlabel('energy error (fraction)'); ylabel('number'); legend(str(9,1:ni));

%axis(h,'tight')
%fcnfontsize(10)
end