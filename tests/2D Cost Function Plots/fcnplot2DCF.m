% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnplot2DCF(input,output,handles,photons,G1)
if output(1).Nsum<19 || output(1).Nsum>1E4;  return;  end;  closeallexcept(handles.GUI.figure1);

xv = linspace(-1,1,100)*input.cube.Lr(1)*.99;
yv = linspace(-1,1,100)*input.cube.Lr(2)*.99;
[X,Y,Z] = meshgrid(xv,yv,G1.p1(1,3));  n=numel(X);

k = fcnoptimizerkr(input,output(1));  k.timeflag=false;  k.reflections=0;
fx = fcnfermatpointvectorized(k,output(1).t,[X(:) Y(:) Z(:) zeros(n,1) ones(n,1)]);  fx=reshape(fx,size(X));

%popoutsubplot(handles.GUI.axes1)
fig; fcnPlotDetector(input,ones(input.cube.pixels,1));
surf(X,Y,Z,fx); shading flat
fcnplot3(G1.p1(1,:),'g+','markersize',15,'linewidth',2); %true point
[~,i]=min3(fx);  fcnplot3([X(i(1),i(2)) Y(i(1),i(2)) G1.p1(1,3)],'ro','markersize',15,'linewidth',2); %min cost
set(gca,'clim',minmax3(fx));






