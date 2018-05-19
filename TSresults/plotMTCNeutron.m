function [] = plotMTCNeutron()
clc; close(findobj(0,'type','figure'));
file = uigetfile();  load(file);  
try
    fname = input.MTC.A.filename;
catch
    fname = MC.FileName;
end

ri=1; %run index
X=MC.xhat(:,:,ri);  T=MC.xtrue(:,:,ri);  MC.elapsedHours
% [r,dx]=fcnrange(X(:,1:3)-X(:,6:8));  dt=X(:,9)-X(:,4);
% ci=~any(isnan(X(:,1:10)),2) & ~all(X==0,2) ...
%     & r>30 ...
%     & dt>2 ...
%     & r./dt>3 ...
%     & X(:,14)>1 ...
% ; 

C = fcngetcandidates(input,[],X,T,[],'MTCNeutron');  ci=all(C,2); %NEUTRON CANDIDATES
%ci=~all(isnan(X),2) & any(X(:,1:3),2); %SIMPLE CANDIDATES
%fprintf('\n     1x     2y     3z     4t     5E     6x     7y     8z     9t    10E    11N    12N  13angle  14np   15v    16aerror 17E0  18E0unc 19Trig  20-24Var(x)\n');
nc=sum(ci);  ne=numel(ci); X=X(ci,:);  T=T(ci,:); [r,dx]=fcnrange(X(:,1:3)-X(:,6:8));  dt=X(:,9)-X(:,4);  E=X-T;

%ha=fig(1,3,1.25); sx={'X','Y','Z'}; for i=1:3; sca; fhistogram(X(:,i),-67:1:67,'FaceColor',fcndefaultcolors(i)); xlabel([sx{i} ' (mm)']); end; title(ha(2),str_(fname));  fcntight;

ha=fig(2,2,1.25);
sca; for i=1:3; fhistogram(X(:,i),-67:2:67); end;   xlabel('Position (mm)');  
sca; for i=1:3; fhistogram(X(:,i)-T(:,i),-67:2:67); end;  xlabel('Error (mm)'); fcntight; legend show

good=false(1536,1); good(fcnpruninglist)=true;
sca; fcnPlotDetector(input,ones(1536,1).*good*.1); fcnplot3(X(:,1:3),'.','Color',fcndefaultcolors(1))
sca;  nb=100;  zv=zeros(nb,nb);  v=linspace(-66,66,nb);
fcnPlotDetector(input,ones(1536,1).*good*.1);  axis off; cla
box=input.cube.box;  handles.detectoroutline = plot3(box.x(:), box.y(:), box.z(:),'Color',[.7 .7 .7]);
[N,xe,ye] = histcounts2(X(:,1),X(:,2),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,ye,zv-67,N,'edgecolor','none');
[N,xe,ye] = histcounts2(X(:,2),X(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(zv+67,xe,ye,N,'edgecolor','none');
[N,xe,ye] = histcounts2(X(:,1),X(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,zv-67,ye,N,'edgecolor','none');  ha(3).CameraViewAngle=8;
fcntight(gca,'csigma'); fcnmarkersize(10); title(ha(1),[repmat(' ',1,80) str_(fname)]); ha(4).CameraViewAngle=8;

%fig; a=double(good); fcnplotdetectorprojection(input,ones(1536,1),a,[0 0 0],'winkeltripel','PMT');  alpha(0.4)

ha=fig(2,3); nb=50;  xe=linspace(0,10,nb);
sca; x=X(:,14); fhistogram(x,.25:.5:5.75); set(gca,'XTick',0:6); xyzlabel('Number of Recoils')
     x=ones(sum(MC.xhat(:,1)==0),1); h=fhistogram(x,.25:.5:5.75); h.FaceColor=[1 1 1]*.7;
sca; x=fcnsigmarejection(X(:,10),6,3); h=fhistogram(x,xe); h.DisplayName=['dE_2' h.DisplayName];
     x=fcnsigmarejection(X(:,5),6,3);  h=fhistogram(x,xe); h.DisplayName=['dE_1' h.DisplayName];
     x=fcnsigmarejection(X(:,17),6,3); h=fhistogram(x,xe); h.DisplayName=['E_0' h.DisplayName]; 
     xyzlabel('E (MeV)','','',sprintf('%g of %g (%.2f) candidate events from %s',nc,ne,nc/ne,str_(fname)));
sca; x=fcnsigmarejection(fcnrange(X(:,1:3)-X(:,6:8)),6,3); fhistogram(x,nb); xyzlabel('dx (mm)');
sca; x=fcnsigmarejection(dt,6,3); fhistogram(x,nb); xyzlabel('dt (ns)');
sca; x=fcnsigmarejection(r./dt,6,3); fhistogram(x,nb); xyzlabel('speed (mm/ns)');
sca; x=(dx./r)*fcnvec2uvec([1 0 0])'; fhistogram(x,linspace(-1,1,nb)); xyzlabel('cos\theta');
%sca; x=fcnsigmarejection(fcnaz(dx)); fhistogram(x,nb);  xyzlabel('azimuth (rad)');
for i=1:6; sca(i); legend show; legend boxoff; end; fcnfontsize(14); fcntight(ha(1:5),'x0')

figure; h=polarhistogram(fcnaz(dx),30); title(str_(fname));

load('neutronconefits.mat')
% F = fits{1,5}; cv=zeros(ne,1); dx=-dx;  %vectors
% %F = fits{2,5}; cv=zeros(ne,1); %fixed cones
% %F = fits{3,5}; cv=real(X(:,13)); %variable cones
%fig(1,3); x=linspace(-1,1,1000); for i=1:3; sca; plot(x,fits{i,5}(x)); end
% for i=[1 3 10]
%     ha=fig(1,2,2);  plotSkyMap(F,dx,cv,i,'neutron');  title(ha(1),str_(input.MTC.A.filename),'Fontsize',16)
% end
for i=[nc]
    ha=fig(1,1,1.5,3);  
    plotSkyMap(fits{1,5},-dx,zeros(nc,1),i,'neutron');  title(ha(1),sprintf('%g of %g events from %s',nc,ne,str_(fname)),'Fontsize',16);
    %plotSkyMap(fits{2,5}, dx,zeros(ne,1),i,'neutron');
    %plotSkyMap(fits{3,5}, dx,real(X(:,13)),i,'neutron');
    ha(1).CameraViewAngle=6;
end
fcnfontsize(12);
export_fig('-a1','-q90','-r300',[fname '.bs.jpg'])


% load M2.mat; input=eval('input');
% load MTCoffsets.mat;  [~,good] = fcnpruninglist;
% fig; fcnPlotDetector(input,ones(1536,1).*good*.1);  view(50,25); xyzlabel('X','Y','Z'); axis on
P1 = mean(X(:,1:3));  %plot3(P1(1),P1(2),P1(3),'.','Color',fcndefaultcolors(1));
P2 = mean(X(:,6:8));  %plot3(P2(1),P2(2),P2(3),'.','Color',fcndefaultcolors(2));  fcnmarkersize(40)
% d = fcnvec2uvec(P2-P1);
% plot3(P2(1)-[0 d(1)*100],P2(2)-[0 d(2)*100],P2(3)-[0 d(3)*100],'-','Linewidth',4)


sprintf('vector = [%.1f %.1f %.1f] +/- [%.1f %.1f %.1f] mm',P2-P1-[2.1 2.9 -2.3],std(X(:,6:8)-X(:,1:3)))
