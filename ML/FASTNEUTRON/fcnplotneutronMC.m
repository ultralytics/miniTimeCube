function [] = fcnplotneutronMC(input,MC,MeVVector)
close all; clc
input.neutron = load('dEfit.mat');
MTCflag = 0;%ischecked(handles.GUI.realdataflag);

h=fig(5,2);
c = {[0 0 1],[0 1 0],[1 0 0],[.8 0 1],[1 .8 0],[0 1 1]};
nb = 50;
x1 = linspace(-1+1/nb,1-1/nb,nb);
av = .5; %alpha value for transparency
ni = numel(MeVVector);
%results.xhat = [xyz t e nphotons, xyz t e nphotons, anglehat, nphat];

name = cell(1,ni);
for i=1:ni
    name{i} = sprintf('%.0fMeV',MeVVector(i));
end
name{1} = 'Bk';
name{2} = 'Pu';

pstr = 'line';
str = cell(8,ni);
fits = cell(3,ni);
for i=1:ni
    if i>6;  c{i}=rand(1,3);  end
    ce = .7*(1-c{i})+c{i}; %edgecolor
    a = zeros(size(MC.xhat(:,1,1),1),13);
    a(:,1) = MC.xhat(:,1,i)~=0 & all(isfinite(MC.xhat(:,:,i)),2); %&  MC.xtrue(:,18,i)>0;
    a(:,2) = MC.xhat(:,6,i)>10; %P1 > 10 photons
    a(:,3) = MC.xhat(:,12,i)>5; %P2 > 2  photons
    a(:,4) = MC.xhat(:,10,i) - MC.xhat(:,4,i) > 1; %dt > 1ns
    a(:,5) = fcnrange(MC.xhat(:,1:3,i) - MC.xhat(:,7:9,i)) > 10; %dx > 10mm
    a(:,6) = MC.xhat(:,5,i)>0.200; %.3 %first bounce > 100keV
    a(:,7) = MC.xhat(:,6,i)./MC.xhat(:,12,i) > .2; %dE1/dE2 must be greather than 0.20
    a(:,8) = MC.xhat(:,5,i)./MC.xhat(:,18,i) < .9; %dE1./E0 < 0.9
    a(:,9) = MC.xhat(:,5,i)./MC.xhat(:,18,i) > .1; %dE1./E0 > 0.1
    a(:,10) = MC.xhat(:,18,i) < 6; %(MeV) E0<20MeV

    %a(:,10) = fcninsidedetector(MC.xhat(:,1:3,i), input, 10); %P1 wall cut
    %a(:,11) = fcninsidedetector(MC.xhat(:,7:9,i), input, 10); %P1 wall cut
    %a(:,12) = MC.collectedPhotonCount(:,2,i)>2; %delayed bookend
    %a(:,13) = all(MC.protonprotonflag(:,:,i),2); %proton-proton bounces only!
    j = all(a(:,1:10),2);   mean(j)
    nj = sum(j);
    mean(a)
    
    continue
    
    %     load net.mat
    %     j = a(:,1)==1;
    %     y = net(MC.xhat(j,:,i)');
    %     j(j) = y>20 & y<70;
    
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
    v0true = -MC.vertex(j,4:6,i);
    
    a=MC.xhat(j,:,i) - MC.xtrue(j,:,i);  if MTCflag; a=MC.xhat(j,:,i); end
    %results.xhat  = [p1 t1 e1 nphotons1, 7p2 t2 e2 nphotons2, 13anglehat, nphat, neutronCapturexyzhat, 18E0hat, E0hatMLpoint, E1hat ];
    %results.xtrue = [p1 t1 e1 nphotons1, 7p2 t2 e2 nphotons2, 13angle,    np   , neutronCapturexyz     18E0     E0          , E1    ];
    
    if nj<10
        for si=1:size(str,1); str{si,i}=''; end
        fprintf('\nWARNING: LESS THAN 10 CANDIDATES FOUND. NOT PLOTTING RESULTS.\n');
        continue
    end
    
    if ~MTCflag
         sca(h(1))
        b=MC.minangleerror(j,1,i); %angles
        y=fcnhist(cos(b*d2r),x1,c{i},'line');
        str{1,2*i-1} = sprintf('%s: %.2f efficiency, %.1f^o 1\\sigma',name{i},sum(j)/numel(j),fcnstdnonlin(b,[],0));
        
        dx = 2/nb; ynormalized = y/sum(y)/dx/(2*pi);  f=createfitexp2(x1, ynormalized);  fits{1,i}=f; %integrate to 1/(2pi)
        plot(x1,f(x1)*dx*nj*2*pi,'-','color',ce,'linewidth',1);
        str{1,2*i} = sprintf('f(x)=%.0gexp^{%.1fx} + %.2fexp^{%.1fx}',f.a,f.b,f.c,f.d);
        
%         %NEW START
%         x1d = linspace(0,pi,50);
%         b=MC.minangleerror(j,1,i); %angles
%         y=fcnhistc(b*d2r,x1d,c{i},pstr);
%         
%         str{1,2*i-1} = sprintf('%s: %.2f efficiency, %.1f^o 1\\sigma',name{i},sum(j)/numel(j),fcnstdnonlin(b,[],0));
%         dx = 2/nb; ynormalized = y/sum(y)/dx/(2*pi);  f=createfitexp2(x1d, ynormalized);  fits{1,i}=f; %integrate to 1/(2pi)
%         plot(x1d,f(x1d)*dx*nj*2*pi,'-','color',ce,'linewidth',1);
%         str{1,2*i} = sprintf('f(x)=%.0gexp^{%.1fx} + %.2fexp^{%.1fx}',f.a,f.b,f.c,f.d);
%         %NEW END
        
        sca(h(2));
        b=a(:,5)./mean(MC.xtrue(j,5,i)); d=fcnsigmarejection(b,6,3); if i==1; x4=linspace(min3(d),max3(d),nb); end
        fcnhist(d,x4,c{i},pstr);
        str{2,i} = sprintf('%s: %.2f 1\\sigma',name{i},fcnstdnonlin(b,[],0));
        
        sca(h(5))
        b=a(:,10); d=fcnsigmarejection(b,3,3); if i==1; x3=linspace(min3(d),max3(d),nb); end
        fcnhist(d,x3,c{i},pstr);
        str{5,i} = sprintf('%s: %.2fns 1\\sigma',name{i},fcnstdnonlin(b,[],0));
        
        sca(h(6))
        b=a(:,7:9); d=fcnsigmarejection(b(:),3,3); if i==1; x2=linspace(min(d),max(d),nb); end
        fcnhist(d,x2,c{i},pstr);
        str{6,i} = sprintf('%s: %.1fmm 1\\sigma',name{i},fcnstdnonlin(b,[],0));
        
        sca(h(3))
        b=a(:,4); d=fcnsigmarejection(b,3,3); %if i==1; x3=linspace(min3(d),max3(d),nb); end
        fcnhist(d,x3,c{i},pstr);
        str{3,i} = sprintf('%s: %.2fns 1\\sigma',name{i},fcnstdnonlin(b,[],0));
        
        sca(h(4))
        b=a(:,1:3); d=fcnsigmarejection(b(:),3,3); %if i==1; x2=linspace(min(d),max(d),nb); end
        fcnhist(d,x2,c{i},pstr);
        str{4,i} = sprintf('%s: %.1fmm 1\\sigma',name{i},fcnstdnonlin(b,[],0));
        
        sca(h(7))
        vt = MC.vertex(j,4:6,i);
        [theta, ct] = fcnangle(p2hat-p1hat,vt);
        y=fcnhistc(ct,x1,c{i},pstr);
        str{7,2*i-1} = sprintf('%s: %.1f^o 1\\sigma',name{i},fcnstdnonlin(theta*180/pi,[],0));
        
        dx = 2/nb; ynormalized = y/sum(y)/dx/(2*pi);  f=createfitrat22(x1, ynormalized);  fits{2,i}=f;
        plot(x1,f(x1)*dx*nj*2*pi,'-','color',ce,'linewidth',1);
        str{7,2*i} = sprintf('f(x)=(%.0gx^2+%.0gx+%.0g)/(x^2+%.0gx+%.0g)',f.p1,f.p2,f.p3,f.q1,f.q2);
        
        sca(h(8)) %true half-angles
        k = p2(:,1)~=0;
        %[theta, ct] = fcnangle(p2(k,:)-p1(k,:),vt(k,:));
        [theta, ct] = fcnangle(p2(k,:)-p1(k,:),p2hat(k,:)-p1hat(k,:));
        y=fcnhistc(ct,x1,c{i},pstr);
        str{8,2*i-1} = sprintf('%s: %.1f^o 1\\sigma',name{i},fcnstdnonlin(theta*r2d,[],0));
        
        dx = 2/nb; ynormalized = y/sum(y)/dx/(2*pi);  f=createfitexp2(x1, ynormalized);  fits{3,i}=f; %integrate to 1/(2pi)
        plot(x1,f(x1)*dx*nj*2*pi,'-','color',ce,'linewidth',1);
        str{8,2*i} = sprintf('f(x)=%.0gexp^{%.1fx} + %.2fexp^{%.1fx}',f.a,f.b,f.c,f.d);
        
        %         sca(h(8))
        %         theta = MC.angleerror(j,:,i);  theta = theta(:);  ct = cosd(theta);
        %         y=fcnhist(ct,x1,c{i},'line');
        %         str{8,2*i-1} = sprintf('%s: %.1f^o 1\\sigma',name{i},fcnstdnonlin(theta,[],0));
        %
        %          dx = 2/nb; ynormalized = y/sum(y)/dx;  f=createfitexp2(x1, ynormalized);  fits{3,i}=f;
        %          plot(x1,f(x1)*dx*nj*361,'-','color',ce,'linewidth',1);
        %          str{8,2*i} = sprintf('f(x)=%.0gexp^{%.1fx} + %.2fexp^{%.1fx}',f.a,f.b,f.c,f.d);
        
        sca(h(9))
        b=a(:,18)./MC.xtrue(j,18,i); d=fcnsigmarejection(b,2,3); if i==1; x7=linspace(min3(d),max3(d),nb); end
        fcnhist(d,x7,c{i},pstr);
        str{9,i} = sprintf('%s: %.2f 1\\sigma',name{i},fcnstdnonlin(b,[],median(b)));
        
        sca(h(10))
        b=a(:,19)./MC.xtrue(j,19,i); d=fcnsigmarejection(b,3,3); if i==1; x8=linspace(min3(d),max3(d),nb); end
        fcnhist(d,x8,c{i},pstr);
        str{10,i} = sprintf('%s: %.2f 1\\sigma',name{i},fcnstdnonlin(b,[],median(b)));
    end
    
    if MTCflag
        load neutronconefits.mat
        v0true = repmat([1 0 0],[nj 1]);
    end
    v21 = p1hat-p2hat;
    skymap(v21,v0true,anglehat,fits{1,i},fits)
end

if ~MTCflag
    sca(h(1))
    xlabel('\Theta'); ylabel('number'); title(sprintf('Neutron Angle Errors\n (%.1fL %s detector, %.1fns, %.2fQE, %.0fyield, %.0fpixels)', ...
        input.cube.volume*1E3,input.Material(1).name,input.scintillatorDecay,input.cube.QEmean,input.Material(1).yield,input.cube.pixels));
    legend(str(1,:));
    sca(h(2))
    xyzlabel('energy error (fraction)','n','','Neutron First Bounce {\bfEnergy Resolution}'); legend(str(2,1:ni));
    sca(h(3))
    title('Neutron First Bounce {\bfTime Resolution}'); xlabel('time error (ns)'); ylabel('number'); legend(str(3,1:ni));
    sca(h(4))
    title('Neutron First Bounce {\bfPosition Resolution}'); xlabel('position error (mm)'); ylabel('number'); legend(str(4,1:ni));
    sca(h(5))
    title('Neutron Second Bounce {\bfTime Resolution}'); xlabel('time error (ns)'); ylabel('number'); legend(str(5,1:ni));
    sca(h(6))
    title('Neutron Second Bounce {\bfPosition Resolution}'); xlabel('position error (mm)'); ylabel('number'); legend(str(6,1:ni));
    
    sca(h(7))
    xlabel('cos(\Theta)'); ylabel('number'); title(sprintf('Estimated P_2 - P_1 vs. True neutron direction')); legend(str(7,:),'Location','best');
    sca(h(8))
    xlabel('cos(\Theta)'); ylabel('number'); title(sprintf('Estimated P_2 - P_1 vs. True P_2 - P_1')); legend(str(8,:),'Location','NorthWest');
    sca(h(9))
    title('E0 = dE_0^1 + E1 {\bfEnergy Resolution}'); xlabel('energy error (fraction)'); ylabel('number'); legend(str(9,1:ni));
    sca(h(10))
    title('E0 = Sum dE all bounces {\bfEnergy Resolution}'); xlabel('energy error (fraction)'); ylabel('number'); legend(str(10,1:ni));
    
    for i=1:numel(h)
        legend(h(i),'boxoff')
    end
    
    axis(h,'tight')
    fcnmarkersize(15)
    fcnlinewidth(1.5)
    %grid(h,'on')
    fcnfontsize(10)
    %x = fcnSNR2angle1sigma(1.8,1)
    save neutronconefits.mat fits
    fprintf('Saving ''neutronconefits.mat'' in ''%s''\n',cd)
end

end

function skymap(v,v0t,angle,f,fits)
nr = 100;
nc = 200;
r = linspace(90,-90,nr);  dr = (r(1) - r(2))*d2r;
c = linspace(-180,180,nc);  dc = (c(2) - c(1))*d2r;
[rm,cm] = meshgrid(r,c);
cc = fcnSC2CCd(1,rm(:),cm(:));  cel=cos(rm(:)'*d2r);
angle = angle*d2r;

h=fig(2,4,1.6,1.2);
nai=0;
for na=[1 3 30 numel(angle)]
    nai=nai+1;
    
    cct = cc';
    v = fcnvec2uvec(v);  v1 = zeros(size(v));
    %na = numel(angle);  %na=30;
    f = fits{2,1};
    t = zeros(na,nr*nc); %theta
    for i=1:na
        if i>1350 
            C1 = fcnVEC2DCM_W2B(isovecs(1));  %background neutron
        else
            C1 = fcnVEC2DCM_W2B(v0t(i,:));
        end
        v1(i,:) = v(i,:)*C1';
        
        ct1 = v1(i,:)*cct;
        %t(i,:) = cos(abs(acos(ct1) - angle(i))); %fits1 method THETA
        t(i,:) = ct1; %fits1 method COSINE THETA
    end
    zm1 = reshape(f(t),size(t));
    
    %NORMALIZE
    nk = 1./( sum(zm1.*cel, 2)' *(dr*dc));

    zmv = nk*zm1; %add
    %zmv = sum(log(nk(ones(nr*nc,1),:)'.*zm1),1);  zmv=exp(zmv-max(zmv));%multiply
    
    
    %PLOT ---------------------------------------------------------------------
    %h=fig(1,2,1.6,1.2);
    %sca(h(1))
    sca(h(nai))
    zm = reshape(zmv,nc,nr);
    pcolor(cm,rm,zm); plot(0,0,'+','MarkerSize',30,'color',[.7 .7 .7]*0);
    shading flat; %ch=colorplot('East'); set(ch,'YColor','k');
    axis equal tight
    xyzlabel('Azimuth (deg)','Elevation (deg)','',sprintf('Neutron Direction Skymap (%.0f neutrons)',i))
    set(gca,'YTick',-90:45:90,'XTick',-180:45:180,'clim',fcnminmax(zm).*[0 1]);
    set(gca,'CameraViewAngle',7,'clim',fcnminmax(zm).*[0 1])
    axis off
    
    sca(h(nai+4))
    fcnPlateCaree2sphere(rm,cm,zm)
    set(gca,'CameraViewAngle',6,'clim',fcnminmax(zm).*[0 1])
    colorbar
end



end





function plotfits
load neutronconefits.mat
fig(1,3,2);
x = linspace(-1,1,1000);
a = acos(x);

sca;  f=fits{1,1};  plot(x,f(x)); title('1')
sca;  f=fits{2,1};  plot(x,f(x)); title('2')
sca;  f=fits{3,1};  plot(x,f(x)); title('3')


end



function fcnPlateCaree2sphere(rm,cm,zm)
sr = 1;
np1 = 600;
np2 = np1+1;
[sx,sy,sz] = sphere(np1);
sc = fcnCC2SC(sx(:),sy(:),sz(:));
el2 = sc(:,2); az2 = sc(:,3);
scL = interp2(rm,cm,zm,reshape(el2,np2,np2),reshape(az2,np2,np2),'linear');
scL = reshape(scL,np2,np2);
sx=sx*sr;
sy=sy*sr;
sz=sz*sr;
surf(sx,sy,sz,scL); shading flat; caxis([min3(zm) max3(zm)+eps])
set(gca,'YDir','Reverse','ZDir','Reverse');
fcn3label; axis off; view(37,30)

fcnplotaxes3(sr*1.25)
fcnplotspherecross(0,0,30,[0 0 0],sr*1.02)
grid off; box on; axis tight equal vis3d;
end

function [fitresult, gof] = createfitexp2(x,y)
[xData, yData] = prepareCurveData(x,y);
ft = fittype( 'exp2' ); % Set up fittype and options.
opts = fitoptions( ft );
opts.Display = 'Off';
opts.StartPoint = [0.687150415982253 1.06648614575917 -0.646882201309983 0.35267553827776];
opts.Normalize = 'on';
[fitresult, gof] = fit( xData, yData, ft, opts ); % Fit model to data.
end


function [fitresult, gof] = createfitrat22(x,y)
[xData, yData] = prepareCurveData(x,y);
ft = fittype( 'rat22' );% Set up fittype and options.
opts = fitoptions( ft );
opts.Display = 'Off';
opts.StartPoint = [0.398589496735843 0.133931250987971 0.0308895487449515 0.939141706069548 0.301306064586392];
opts.Normalize = 'on';
[fitresult, gof] = fit( xData, yData, ft, opts ); % Fit model to data.
end
