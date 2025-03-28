% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

clc; close; format long g

MeV = 11;
GEANT = fcnloadGEANT(['MTCG4SimOutput_ReactorPositronAndNeutronInEJ254WithNaturalBoronDopedAt5PercentLargeCube_' num2str(MeV) 'MeV_10000runs.dat']);

zv1     = zeros(1E4,1);
zv3     = zeros(1E4,3);
ppN     = zv3;
np1     = zv3;
puVec   = zv3;
nuVec   = zv3;
pE      = zv1;
nE      = zv1;

rpy     = [0 90 0]*pi/180;

e=0;
for i = 1:GEANT.ne
    input.eventNumber=GEANT.ve(i);
    G1 = fcnGEANT2G1(GEANT,input);
    v1=find(G1.tid==1);  v2=find(G1.tid==2);
    
    if numel(v1)>2 && numel(v2)>2
        e=e+1;

        ppN(e,1:3)      = G1.p1(v1(end),:); %last positron point
        np1(e,1:3)      = G1.p1(v2(2),:); %first neutron bounce
        
        puVec(e,1:3)    = G1.uvec.tid1;
        nuVec(e,1:3)    = G1.uvec.tid2;
        pE(e,1)         = max(G1.ke1(v1));
        nE(e,1)         = max(G1.ke1(v2));
    end
    fprintf('event %.0f\n',i)
end
ne=e; %number of events with at least two points in positron streak
ppN=ppN(1:e,:); np1=np1(1:e,:); puVec=puVec(1:e,:); nuVec=nuVec(1:e,:); pE=pE(1:e); nE=nE(1:e);

% E = sqrt(P^2*c^2+(m*c^2)^2)
% P = sqrt(E^2+(m*c^2)^2)/c
% E1^2 - (m1*c^2)^2  =  E2^2 - (m2*c^2)^2

mass.neutron = 939.565346; %MeV/c^2
mass.positron = .510998918; %MeV/c^2
mass.proton = 938.272013; %MeV/c2
c = 299792458; %m/s

constant    = 0.782588081999961; % mass.neutron-mass.positron-mass.proton

mode        = 4;
v           = zeros(ne,3);

nx          = 100;
ny          = 100;
sx          = linspace(0,5,nx); %endpoint position sigma
sy          = linspace(0,100,ny); %photon count sigma
SNR1       = zeros(nx,ny);
SNR2     = zeros(nx,ny);

SNR1         = zeros(nx,ny);
SNR2         = zeros(nx,ny);
p1          = zeros(ne,3);
p2          = zeros(ne,3);
for i = 1:nx
    for j = 1:ny
        x1 =       randn(ne,3)*sx(i);  %[0 0 0]
        x2 = ppN + randn(ne,3)*sx(i);
        x3 = np1 + randn(ne,3)*sx(i);
        
        puvec2 = fcnvec2uvec(x2-x1);
        nuvec2 = fcnvec2uvec(x3-x1);
        
        pE2 = sqrt(max((pE+mass.positron + randn(ne,1)*sy(j)/input.fluid.yield).^2-mass.positron^2,1/input.fluid.yield));
        nE2 = sqrt(max((nE+mass.neutron + randn(ne,1)*sy(j)/input.fluid.yield).^2 + constant^2 - mass.neutron^2,1/input.fluid.yield));
        for k=1:3
            p1(:,k) = pE2.*puvec2(:,k);
            p2(:,k) = nE2.*nuvec2(:,k);
        end

        v = fcnvec2uvec(p1+p2);  
        s=mean(v);  s=s(1);  n=mean(std(v));  SNR1(i,j)=s/n;
       
        
        v = fcnvec2uvec(x3-x1);  
        s=mean(v);  s=s(1);  n=mean(std(v));  SNR2(i,j)=s/n;
    end
end

p1str = 'Point 1  =  antineutrino vertex';
p2str = 'Point 2  =  positron annihilation point';
p3str = 'Point 3  =  first neutron bounce';

% %MOMENTUM RESULTS ---------------------------------------------------------
% set(figure,'Position',[10 40 500 500]); z = azstd1';
% pcolor(sx,sy,z); xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); shading flat
% title({[sprintf('%.0fMeV Momentum Unit Vectors (%.0f)',MeV,ne) ' Azimuth Angle 1\sigma'], p1str, p2str, p3str})
% colormap(hot); colorbar; caxis([0 122]); hold on
% [~, h] = contour3(sx,sy,SNR1'-SNR2',[0 0],'g'); set(h,'LineWidth',2)
% legend('Az Std (deg)','Advantage Boundary','CHOOZ Az Std = 119.65^o')
% % Create the textarrow object: 
% x=[.2 .14];  y=[.15 .12];  annotation('textarrow',x,y,'string',['Min=' sprintf('%.2f',min(min(z))) '^o'],'FontSize',14,'Color',[.7 .7 .7]);
% 
% %RECONSTRUCTION VECTOR RESULTS --------------------------------------------
% set(figure,'Position',[120 40 500 500]); z = azstd2';
% pcolor(sx,sy,z); shading flat; xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); 
% title({[sprintf('%.0fMeV Momentum Unit Vectors (%.0f)',MeV,ne) ' Azimuth Angle 1\sigma'], p1str, p3str})
% colormap(hot); colorbar; caxis([0 122]); hold on
% [~, h] = contour3(sx,sy,SNR1'-SNR2',[0 0],'g'); set(h,'LineWidth',2)
% legend('Az Std (deg)','Advantage Boundary','CHOOZ Az Std = 119.65^o')
% % Create the textarrow object: 
% x=[.2 .14];  y=[.15 .12];  annotation('textarrow',x,y,'string',['Min=' sprintf('%.2f',min(min(z))) '^o'],'FontSize',14,'Color',[.7 .7 .7]);
% 
% %BOTH RESULTS OVERLAID ----------------------------------------------------
% set(figure,'Position',[230 40 900 700])
% h = mesh(sx,sy,SNR1',ones(nx)*20);  hidden off; xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); zlabel('Azimuth Angle 1\sigma (deg)')
% caxis([0 122]); hold on
% mesh(sx,sy,SNR2',ones(nx)*100); %hidden off
% legend('Momentum Equation','Reconstruction Vector')
% title({[sprintf('%.0fMeV Momentum Unit Vectors (%.0f)',MeV,ne) ' Azimuth Angle 1\sigma'], p1str, [p2str '(Momentum Equation ONLY)'], p3str})
% set(gca,'ZLim',[0 121])
% %plots2PDFs


%MOMENTUM RESULTS ---------------------------------------------------------
zm = max3([SNR1 SNR2]); 
h=fig(2,2);
axes(h(1))
z = SNR1';
pcolor(sx,sy,z); xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); shading flat
title(sprintf('%.0fMeV Momentum Equation Vector SNR',MeV))
colorbar; caxis([0 zm]); hold on
[~, hc] = contour3(sx,sy,SNR1'-SNR2',[0 0],'g'); set(hc,'LineWidth',2)
legend('vector SNR','Advantage Boundary')
x=[.1 .06];  y=[.57 .55];  annotation('textarrow',x,y,'string',['max=' sprintf('%.2f',max3(z))],'FontSize',14,'Color',[.7 .7 .7]);

%RECONSTRUCTION VECTOR RESULTS --------------------------------------------
axes(h(2))
z = SNR2';
pcolor(sx,sy,z); xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); shading flat
title(sprintf('%.0fMeV CE2CE Vector SNR',MeV))
colorbar; caxis([0 zm]); hold on
[~, hc] = contour3(sx,sy,SNR1'-SNR2',[0 0],'g'); set(hc,'LineWidth',2)
legend('vector SNR','Advantage Boundary')
x=[.6 .56];  y=[.57 .55];  annotation('textarrow',x,y,'string',['max=' sprintf('%.2f',max3(z))],'FontSize',14,'Color',[.7 .7 .7]);
fcnfontsize(8)

%BOTH RESULTS OVERLAID ----------------------------------------------------
axes(h(3))
mesh(sx,sy,SNR1',ones(nx)*20); hold on; hidden off; xlabel('Position Estimate 1\sigma (mm)'); ylabel('Photon Count Estimate 1\sigma (#)'); zlabel('SNR')
mesh(sx,sy,SNR2',ones(nx)*100); caxis([0 120]); set(gca,'zlim',[0 zm])%hidden off
legend('Momentum Equation','Reconstruction Vector','location','northeast')
title('Overlay Comparison')
view(25,35)
%plots2PDFs

delete(h(4))
