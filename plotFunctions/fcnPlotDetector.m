% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function h = fcnPlotDetector(input,Cdata,AlphaData)
n = input.cube.pixels;

if nargin < 2
    Cdata = zeros(n,1);
end

if nargin < 3
    AlphaData = ones(n,1)*.08;
    i = Cdata>max(Cdata)*.01;
    AlphaData(i) = 0.8*.5;
    Cdata(~i) = 0;
end

try
    flags = evalin('base','flags');
    af = find(flags.status.plotside); %active faces
catch
   af=[1 2 3 4 5 6]; 
end

i = find(any(input.cube.all.fi==af,2)); %active pixel indices

h=patch(input.cube.all.v.x(i,:)', input.cube.all.v.y(i,:)', input.cube.all.v.z(i,:)','w', ...
    'EdgeColor','none', ... 'EdgeAlpha',0, ...
    'FaceColor','Flat', ...
    'FaceAlpha','Flat', ...
    'FaceVertexCData',Cdata(i), ...
    'FaceVertexAlphaData',AlphaData(i), ...
    'CDataMapping','scaled', ...
    'AlphaDataMapping','none', ...
    'FaceLighting','none', ...
    'BackFaceLighting','unlit', ...
    'parent',gca);

axis tight equal off
set(gca,'ydir','reverse','zdir','reverse','CameraViewAngle',9)
fcnview('skew')
fcntight('csigma')

%PLOT BOX
box=input.cube.box;  handles.detectoroutline = plot3(box.x(:), box.y(:), box.z(:),'Color',[.7 .7 .7]);

% %PLOT NTC ASICs
% S=NTCIRSDmap(1:128,'SRCCH'); A=S(:,2)*4+S(:,3);
% fig; fcnPlotDetector(input,A,ones(128,1));  h=colorbar('East'); h.Label.String='ASIC';


% %PLOT PIXEL NUMBERS
% xyz=input.cube.all.xyz;
% for i=1:input.cube.pixels
%     text(xyz(i,1),xyz(i,2),xyz(i,3),sprintf('%g',i),'color',[1 1 1]*.3,'fontsize',10,'horizontalalignment','center','verticalalignment','middle');
% end
