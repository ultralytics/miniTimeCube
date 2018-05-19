function [] = plotTSMTC()
clc; close(findobj(0,'type','figure'));

file = uigetfile();
load(file);  A=input.MTC.A;

f=1;
X=MC.xhat(:,:,f); T=MC.xtrue(:,:,f);   
ci=~any(isnan(X(:,1:5)),2) & ~all(X==0,2);  ne=sum(ci);  X=X(ci,:);  T=T(ci,:);
MC.collectedPhotonCount=MC.collectedPhotonCount(ci,:); 
MC.vertex = MC.vertex(ci,:);


%[~,F]=polyfitrobust(X(:,5),T(:,5),5,2.5,6,'horizontal');  save MTCEcalCo60.mat F
load MTCEcalPoint.mat;  X(:,5)=F(X(:,5));  %635 Pixel MTC ENERGY CALIBRATION

%neta = NNtrain(X(:,[1]),T(:,[1]),10,1); % save('NN.MTC Calibration 0-10 MeV.mat','neta')
%X(:,[1])=neta(X(:,[1])')';

%fig; plot(X(:,5),T(:,5),'.','Markersize',4);  x=linspace(0,max(X(:,5)),1000); 
%plot(polyval(Pe,x),x,'-','DisplayName','Polynomial');
%plot(x,neta(x),'-','DisplayName','Neural Network'); fcnlinewidth(2)

%X=T;

E=X-T;  xp=X(:,1:5);  tp=T(:,1:5); if size(X,2)>=10; xd=X(:,(1:5)+5); td=T(:,(1:5)+5);  n=2; else n=1; end; n=1;
load M2.mat; input=eval('input');
load MTCoffsets.mat;  [~,good] = fcnpruninglist;
for i=1:n
    switch i
        case 1;  x=xp;  t=tp;  s='PROMPT';
        case 2;  x=xd;  t=td;  s='DELAYED';
    end
    ha=fig(2,2,1,1.1);  c=fcndefaultcolors(i);  ne=sum(~any(isnan(x),2) & ~all(x==0,2));

    nb=100;
    sca; histogram(x(:,4),nb,'facecolor',c); xyzlabel('T (ns)');  title(sprintf('%g %s Events\n%s',ne,s,MC.FileName)); ha(3).Title.FontSize=16;
    if any(regexpi(MC.particleName{:},'neutrino'))
        sca; a=x(:,5); [hy,hx]=histcounts(a,linspace(0,max(fcnsigmarejection(a,3,3)),nb)); xyzlabel('E (MeV)');
        hx=hx(2:end)/2+hx(1:end-1)/2;  g=fcnspecreactoremission(hx+1.78).*fcnspecdetection(hx+1.78);  g=g./max(g);  if i==1; hy=hy.*g; end
        bar(hx,hy,1,'FaceColor',c,'FaceAlpha',.6);
        NS{i}={hx,hy};  save simulatedNuebarSpectrum.mat NS %#ok<AGROW>
    else
        sca; a=x(:,5); h=histogram(a,linspace(0,max(fcnsigmarejection(a,4,3)),nb),'facecolor',c,'DisplayName','Fit'); xyzlabel('E (MeV)');
        %if ~all(t(:)==0); a=t(:,5); histogram(a,h.BinEdges,'DisplayName','True','facecolor',[.7 .7 .7]); xyzlabel('E (MeV)'); end
    end
    
    nb=300; 
    sca;  zv=zeros(nb,nb);  v=linspace(-66,66,nb);
    fcnPlotDetector(input,ones(1536,1).*good*.1);  axis off; cla
    box=input.cube.box;  handles.detectoroutline = plot3(box.x(:), box.y(:), box.z(:),'Color',[.7 .7 .7]);
    [N,xe,ye] = histcounts2(x(:,1),x(:,2),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,ye,zv-67,N,'edgecolor','none');
    [N,xe,ye] = histcounts2(x(:,2),x(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(zv+67,xe,ye,N,'edgecolor','none');
    [N,xe,ye] = histcounts2(x(:,1),x(:,3),v,v);  [xe,ye]=ndgrid(xe,ye);   surf(xe,zv-67,ye,N,'edgecolor','none');
    fcntight('csigma');
    sca; plot3(x(:,1),x(:,2),x(:,3),'.','color',c); fcnPlotDetector(input,ones(1536,1).*good*.1);  axis off
    
    ha(3).CameraViewAngle=7; ha(4).CameraViewAngle=7; %for i=1:3; ha(i).XLim=[-67 67]; ha(i).YLim=[-67 67]; axis(ha(i),'square'); ha(i).YTick=ha(i).XTick;  box(ha(i),'on'); grid(ha(i),'off'); end; axis(ha(4:end),'tight'); fcntight(ha(1:3),'c')
    fcnmarkersize(1.5); fcntight; fcntight('csigma')
end
fig(1,3,1.25); nb=100; a=xp(:,5);  i=a<2;  labels={'E (MeV)','Triggered Channels per event','PEs per event'};
h=sca; histogram(a,linspace(0,max(a(i)),nb),'facecolor',c,'DisplayName','Fit'); xyzlabel(labels{1}); h.YScale='log'; h.YLim(1)=.5;
h=sca; histogram(double(MC.vertex(i,10)),nb); xlabel(labels{2}); h.YScale='log'; h.YLim(1)=.5; title(MC.FileName)
h=sca; histogram(double(MC.collectedPhotonCount(i,1)),nb); xlabel(labels{3}); h.YScale='log'; h.YLim(1)=.5;


% %MTC CF252 SPECTRA
% clear a b y x;
% a={xp(:,5),MC.vertex(:,10),MC.collectedPhotonCount(:,1),MC.collectedPhotonCount(:,2)};
% x={linspace(0,2,nb),linspace(0,400,nb),linspace(0,500,nb),linspace(0,15,nb)};
% for i=1:4
%     y{i}=histcounts(a{i},x{i}); %xyzlabel('E (MeV)'); h.YScale='log'; h.YLim(1)=.5;
% end
% %save overlay2747.mat y
% 
% 
% fig(1,3,1.25);
% b=load('overlay2742');
% for i=1:3
%     sca;
%     z=fcnsmooth(b.y{i}',5);                     plot(edges2centers(x{i}),z/max(z));
%     z=fcnsmooth(y{i}',5)*250/105;               plot(edges2centers(x{i}),z/max(z));
%     z=fcnsmooth(y{i}'*250/105-b.y{i}',5);       plot(edges2centers(x{i}),z/max(z)); 
%     %z=(y{i}'*250/105-b.y{i}')./double(b.y{i}');     plot(edges2centers(x{i}),max(fcnsmooth(z,3),0)); 
%     xlabel(labels{i})
%     fcntight('y0')
% end
% fcnlinewidth(2)
% sca; legend('Background','Cf252','Cf252 - Background','Signal to Background Ratio')
% 
% fig; fhistogram(MC.collectedPhotonCount(:,2),linspace(0,15,200)); xlabel('timestamp \sigma'); title(MC.FileName);


ha=fig(2,4,1);  nb=90;  x=linspace(-67,67,nb);  t='tile';
sca; y=fcnhistc(E(:,1),x,fcndefaultcolors(1)); xyzlabel('X error (mm)');  %histNfit(x,y,'b');
sca; y=fcnhistc(E(:,2),x,fcndefaultcolors(2)); xyzlabel('Y error (mm)');  %histNfit(x,y,'b');
sca; y=fcnhistc(E(:,3),x,fcndefaultcolors(3)); xyzlabel('Z error (mm)');  %histNfit(x,y,'b');
sca; fcnhistc(E(:,5),nb,fcndefaultcolors(4));  legend show; xyzlabel('E error (MeV)');
sca; hist2(X(:,1),T(:,1),{x,x},t); xyzlabel('X fit (mm)','true')
sca; hist2(X(:,2),T(:,2),{x,x},t); xyzlabel('Y fit (mm)','true')
sca; hist2(X(:,3),T(:,3),{x,x},t); xyzlabel('Z fit (mm)','true')
sca; a=[X(:,5) T(:,5)]; x=linspace(min3(a),max3(a),nb); hist2(X(:,5),T(:,5),{x,x},t); xyzlabel('E fit (MeV)','true')
axis(ha,'tight')
end