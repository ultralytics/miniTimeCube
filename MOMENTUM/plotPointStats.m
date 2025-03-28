% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = plotPointStats(G)

p1 = zeros(1,3);
p2 = zeros(1,3);
d0 = zeros(1,3);
d1 = zeros(1,3);
d2 = zeros(1,3);

%COLLECT DATA -------------------------------------------------------------
iv=0;
yield = 1;
BQF = 1;
for i = 1:10000
    if G.positron(i).eIoni.count>1 && G.neutron(i).hIoni.count>1
        iv = iv+1;
        %POSITRON
        uvec = G.positron(i).uVec.CC;
        C = fcnVEC2DCM_W2B(uvec);
        p1(iv,1:3) = (C*G.positron(i).eIoni.CC(1,1:3)')';
        p2(iv,1:3) = (C*G.positron(i).eIoni.CC(2,1:3)')';
        p1(iv,4)   = G.positron(i).eIoni.dE(1)*yield;
        p2(iv,4)   = G.positron(i).eIoni.dE(2)*yield;
        
        %NEUTRON
        uvec = G.neutron(i).uVec.CC;
        C = fcnVEC2DCM_W2B(uvec);
        d0(iv,1:3) = (C*G.neutron(i).p1.CC(1,1:3)')';
        d1(iv,1:3) = (C*G.neutron(i).hIoni.CC(1,1:3)')';
        d2(iv,1:3) = (C*G.neutron(i).hIoni.CC(2,1:3)')';
        
        d0(iv,4)   = 0;
        d1(iv,4)   = G.neutron(i).hIoni.dE(1)*yield*BQF;
        d2(iv,4)   = G.neutron(i).hIoni.dE(2)*yield*BQF;
    end
end


%PLOT RESULTS -------------------------------------------------------------
dataCells = {p1,p2,d0,d1,d2};
c = {'r','r','m','b','b'};
titleCells = {' 1^s^t Positron eIoni Points', ' 2^n^d Positron eIoni Points' ...
    ' 1^s^t Neutron Points (NO dE!!)', ' 1^s^t Neutron hIoni Points', ' 2^n^d Neutron hIoni Points'};
    
alpha = 1;
for ix = 1:5
    x = dataCells{ix};
    figure
    set(gcf,'Position',[10+100*(ix-1) 40 700 700]);
    
    h = subplot(2,2,1); set(h,'ZDir','Reverse','YDir','Reverse')
    minx = min([x(:,1); -2]);
    plot3([minx 0],[0 0],[0 0],'g','LineWidth',3); hold on
    plot3(x(:,1),x(:,2),x(:,3),'.','MarkerSize',1,'Color',c{ix}); view(-150,50)
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); title({[num2str(iv) titleCells{ix}], ['WRT ' num2str(G.neutrino.MeV(1,1)) 'MeV Neutrino Annihilation Point'],['Expressed in ' num2str(iv) ' Respective Truth Vector Frames']})
    axis tight equal vis3d
    box on
    
    subplot(2,2,2);
    hist(x(:,4),30); h = get(gca,'Children'); axis tight
    set(h,'FaceColor',c{ix},'EdgeColor',c{ix},'FaceAlpha',alpha,'EdgeAlpha',alpha);
    xlabel('Proton dE (MeV)'); ylabel('Number of Events'); title({'Proton dE Produced by', titleCells{ix}})
    
    subplot(2,2,3)
    r = normRows(x);
    hist(r,30); h = get(gca,'Children'); axis tight
    set(h,'FaceColor',c{ix},'EdgeColor',c{ix},'FaceAlpha',alpha,'EdgeAlpha',alpha);
    xlabel('Range (mm)'); ylabel('Number of Events'); title({['Ranges of ' titleCells{ix}], ' From Neutrino Annihilation Point'})
    
    subplot(2,2,4)
    theta = thetaRows(x);
    hist(theta,linspace(0,180,30)); h = get(gca,'Children'); axis tight
    set(h,'FaceColor',c{ix},'EdgeColor',c{ix},'FaceAlpha',alpha,'EdgeAlpha',alpha);
    xlabel('Angle (deg)'); ylabel('Number of Events'); title({['Angle Of ' titleCells{ix} ], 'Off Of Truth Vector'})
end

end

function r = normRows(x)
r = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2)+eps;
end

function theta = thetaRows(v2)
n = length(v2);
v1 = ones(n,1)*[1 0 0];
r = normRows(v2);
v2 = [v2(:,1)./(r+eps)  v2(:,2)./(r+eps)  v2(:,3)./(r+eps)];
theta=acos( v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3))*180/pi;
end
