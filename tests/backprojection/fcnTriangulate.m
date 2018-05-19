function [voxels, error, handlesOut] = fcnTriangulate(input,output,handles,photons)
error.flags = [output(1).Nsum<10   output(2).Nsum<3      0]; %last flag is fminsearch 'exitflag'
handlesOut = [];
if output(1).Nsum==0
    voxels = [];
    return
end
tic

%CONSTANTS ----------------------------------------------------------------
zt              = output(1).t; %timestamp measurements
s               = 0.100; %timestamp noise (ns)

%USER INPUTS --------------------------------------------------------------
nv1     = 100; %number of voxels per side
nt      = 240; %number of time frames
tvec    = linspace(-.8,.8,nt);
%tvec    = .05;

%BEGIN --------------------------------------------------------------------
[X,Y,Z]     = ndgrid(linspace(-1,1,nv1)*input.cube.Lr(1)*1, linspace(-1,1,nv1)*input.cube.Lr(1)*1, 0);
ne          = numel(X); %number of elements
%[x, y, z]  = fcnpolyorthostreakfit(input, photons);
%DCM        = fcnVEC2DCM_B2W([x(2)-x(1) y(2)-y(1) z(2)-z(1)]);
coeff = pca(photons.sourcePos);
DCM = fcnVEC2DCM_B2W(coeff(:,1));
voxels.pos = [reshape(X,ne,1) reshape(Y,ne,1) reshape(Z,ne,1)]*DCM' + ones(ne,1)*mean(photons.sourcePos);
X=reshape(voxels.pos(:,1),nv1,nv1); Y=reshape(voxels.pos(:,2),nv1,nv1); Z=reshape(voxels.pos(:,3),nv1,nv1);


%NEW APPROACH
p = voxels.pos;
k = fcnoptimizerkr(input,output(1));

dx = k.pxyz(:,1)-p(:,1)';  dy = k.pxyz(:,2)-p(:,2)';  dz = k.pxyz(:,3)-p(:,3)';
dotprod = k.nxyz(:,1).*dx + k.nxyz(:,2).*dy + k.nxyz(:,3).*dz;
rs      = dx.^2+dy.^2+dz.^2;
r      	= sqrt(rs); %(mm)
ct    	= dotprod./r;  ct=ct(k.pid,:); %cos(theta)

t = -r(k.pid,:)*k.ci + zt;

method = 'faster';
pdfmethod = 'analytical';
switch method
    case 'lessmemory' %uses less ram but 2X slower
        v = cell(nt,1);
        for j=1:nt
            fprintf('%.0f/%.0f\n',j,nt)
            dt = t-tvec(j);
            switch pdfmethod
                case 'LUT'
                    v{j} = sum( interp1c(k.smearExp.x,k.smearExp.ys,dt) );
                case 'analytical'
                    s = .02; %ns
                    %                     pdf0 = 1/(2*pi*s^2) * exp(-dt.^2/(2*s^2));
                    %                     pdf1 = fcnpdfnorm(t,tvec(j),s);
                    %                     pdf2 = fcnpdfEMG(t,tvec(j),s,1/.01);
                    v{j} = (1/(2*pi*s^2)) * sum( exp(-dt.^2/(2*s^2)) );
            end
        end
        v = cell2mat(v)';
    case 'faster' %faster but may use up all ram!
        dt = t + reshape(-tvec,[1 1 nt]);
        switch pdfmethod
            case 'LUT'
                v = sum( nterp1c(k.smearExp.x,k.smearExp.ys,dt) );
                
            case 'analytical'
                %v = sum(fcnpdfnorm(dt,0,s));
                v = (1/(2*pi*s^2)) * sum( exp(-dt.^2/(2*s^2)) );
        end; clear dt
end
v = reshape(v,[nv1 nv1 nt]);

fprintf('Backscatter Complete in %.1fs, saving as ''%s''...',toc,'backscatter.mat');
save -v6 backscatter.mat v tvec nt nv1 X Y Z
fprintf(' Done.\n')

x = linspace(1,nv1,nv1*2);
v = interpn(v,x',x',linspace(1,nt,nt*2),'linear'); nt=nt*2; tvec=linspace(min(tvec),max(tvec),nt)
X=interpn(X,x,x'); Y=interpn(Y,x,x'); Z=interpn(Z,x,x');
%rho = fcndsearch(v, .90);

%MAKE VIDEO!
[h, hf] = fig(1,2,1920,1080);
popoutsubplot(handles.GUI.axes1,h(1));
popoutsubplot(handles.GUI.axes1,h(2));
box(h(1),'on'); box(h(2),'on'); axis(h,'vis3d','off'); h1=[]; a=.75; hold(h(1),'on'); hold(h(2),'on')

V=max(v,[],3);
surf(h(2), X, Y, Z, V,'EdgeColor','none','FaceAlpha',a);  caxis(h(2),[min3(V) max3(V)+eps]);
title('Backprojection Maximum')

hl = findobj(gcf,'type','legend'); set(hl,'location','northeast'); colormap(jet)

fcncolorbar(a,[],hf); caxis(h(1),[min3(v) max3(v)+eps]);

%return
fname=fcnincrementfname('backprojection.mp4'); vidObj=VideoWriter(fname,'MPEG-4'); vidObj.Quality=100; open(vidObj);
%h1=surf(h(1), X, Y, Z, Z*0,'EdgeColor','none','FaceAlpha',a);
for i = 240%1:nt
    vi = v(:,:,i);
    deleteh(h1); %[az,el]=view; view(az+1/nt*20,el-1/nt*10)
    h1=surf(h(1), X, Y, Z, vi,'EdgeColor','none','FaceAlpha',a);
    %set(h1,'Cdata',vi);
    title(h(1),sprintf('Backprojection for t=%.3fns',tvec(i)))
    drawnow('update'); 
    writeVideo(vidObj,getframe(hf));
end
fprintf('Video Complete, saving as ''%s''...',fname); close(vidObj); fprintf(' Done.\n')




