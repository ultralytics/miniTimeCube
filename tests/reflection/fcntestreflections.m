% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcntestreflections(input,output,handles,G1,photons)
clc; close all
ha = fig(1,2,2.5,2); sca(ha(1));
%copyobj(get(handles.GUI.axes1,'children'), ha(1)); 
xyzlabel('X (mm)','Y (mm)','Z (mm)'); set(gca,'YDir','Reverse','ZDir','Reverse'); axis tight equal% vis3d
%[az1, el1] = view(handles.GUI.axes1);  view(ha(1), az1,el1);  set(gca,'CameraViewAngle',7.5)
cube(input.cube.Lr,[0 0 0],1);

p = G1.p1(1,:);
p = [300 200 0];
t = G1.t1(1);
X = [p 0 1];
%X = [X [-40 40 40 0 1] ]; %second point

vidflag = false;
if vidflag;  fname=fcnincrementfname('Reflection Model.mp4'); vidObj=VideoWriter(fname,'MPEG-4');  set(vidObj,'Quality',100);  open(vidObj);  end %#ok<*USENS>
tic

fi = [1 2 3 4 5 6]; %[2 3 6];% goes with skew view; face indices
input.Material(1).mu(1) = 1.58;
input.Material(1).mu(2) = 1.31; 
iap = find(any(input.cube.all.fi==fi,2));
k = fcnoptimizerkr(input,output(1));  k.timeflag = 0;  k.QEmap=1;

sca(ha(1))
hs = patch(input.cube.all.v.x(iap,:)', input.cube.all.v.y(iap,:)', input.cube.all.v.z(iap,:)','w', ...
    'EdgeColor','none', ...
    'EdgeAlpha',0, ...
    'FaceColor','Flat', ...
    'FaceAlpha','Flat', ...
    'FaceVertexCData',zeros(k.npixels,1), ...
    'FaceVertexAlphaData',.7, ...
    'CDataMapping','scaled', ...
    'AlphaDataMapping','none', ...
    'FaceLighting','none', ...
    'BackFaceLighting','unlit', ...
    'Parent',ha(1));  view(-150,26);  colorbar('location','East'); axis off

mmc = [0 eps];
for t = linspace(0,1,1);
    zt = ones(k.npixels,1)*t;
    
    for ii=1:3
    [~, ~, ~, ~, ft] = fcnfermatpointvectorized(k,zt,X);
    end
    set(hs,'FaceVertexCData',ft);
    title(ha(1),sprintf('n_1=%.2f, n_2=%.2f, t=%.3fns',k.ior1,k.ior2,t))
    [~, ~, ~, mmct] = fcnaxesdatalims(ha(1));  mmc(2) = max(mmc(2),mmct(2));  set(ha(1),'clim',mmc);

    sca(ha(2))
    fcnplotdetectorprojection(input,iap,ft,p*0); axis off
    
    %fcntight('c')
    drawnow
    if vidflag; pause(.1);  frame=getscreen(gcf);  writeVideo(vidObj,frame);  end %#ok<*UNRCH>
end
toc
sca(ha(1)); colormap(jet(256))
alpha(.7)
[sum(ft) max(ft) mean(ft)]
if exist('vidObj','var');  close(vidObj); end %#ok<NODEF>
%sca(handles.GUI.axes1)
%c = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 .5 0; 1 0 0; 1 0 1];
%c = fcndefinecolormap(c);
%c = jet;
%colormap(c)

h=fcnplotline(p,'w.-'); set(h,'markersize',30,'linewidth',2,'color',[0 0 1])
% ha = fig(1,2,1.95,2);
% popoutsubplot(hsp,ha(1)); sca(ha(1)); title('NEAR'); axis off; view(-90,0); colormap(c); colorbar; set(gca,'clim',fcnminmax(fxv{3}) + [0 eps])
% popoutsubplot(hsp,ha(2)); sca(ha(2)); title('FAR');  axis off; view( 90,0); colormap(c); colorbar; set(gca,'clim',fcnminmax(fxv{1}) + [0 eps])
% %set(ha,'CameraViewAngle',1.5)

end