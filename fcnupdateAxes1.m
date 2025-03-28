% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [input,output,G1,photons,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles)

if flags.update.edit1 %NEW EVENT
    sca(handles.GUI.axes1);
    [input,output,G1,photons,PE,flags] = fcn1Event(input,flags,handles,GEANT);
    
    %SET UP TRACK ID PLOTTING HANDLES
    deleteh(handles.upid);      handles.upid=gobjects(1,G1.nupid);
    deleteh(handles.upidsp);    handles.upidsp=gobjects(1,G1.nupid);
    deleteh(handles.upidcp);    handles.upidcp=gobjects(1,G1.nupid);
    deleteh(handles.gridanode); 
    deleteh(handles.stripanode);
    deleteh(handles.segments);
    for i=1:G1.nupid
        handles.upid(i)   = plot3(nan,nan,nan,'.-','color',G1.upidcolor{i},'MarkerSize',30,'Display',G1.upidname{i},'LineWidth',2.5); %GEANT particles
        handles.upidsp(i) = plot3(nan,nan,nan,'-','color',G1.upidcolor{i},'Display','','Parent',handles.GUI.axes1,'LineWidth',3); %scintillation photons
        handles.upidcp(i) = plot3(nan,nan,nan,'-','color',G1.upidcolor{i}*.5,'Display','','LineWidth',3); %cherenkov photons
    end
    
    handles.gridanode = patch(1, 1, 1, 'w', ...
        'EdgeColor','none', ...
        'FaceColor','Flat', ...
        'FaceAlpha','Flat', ...
        'FaceVertexCData',1, ...
        'AlphaDataMapping','none', ...
        'FaceVertexAlphaData',1, ...
        'parent',handles.GUI.axes1);
    handles.segments = patch(1, 1, 1, 'w', ...
        'EdgeColor',[.5 .5 1], ...
        'FaceColor',[1 1 1]*.9, ...
        'FaceAlpha',.3, ...
        'parent',handles.GUI.axes1);
    if isfield(photons,'fibervx'); handles.segments.XData=photons.fibervx; handles.segments.YData=photons.fibervy; handles.segments.ZData=photons.fibervz; end
    handles.stripanode = surface(1,1,1,0,'facecolor','flat','edgecolor','none','FaceAlpha',.6,'Parent',handles.GUI.axes1);
    
    
    %DETECTOR OUTLINE
    deleteh([handles.detectoroutline(:);  handles.detectorfacelabels(:)]);
    box=input.cube.box;  handles.detectoroutline = plot3(box.x(:), box.y(:), box.z(:),'Color',[.7 .7 .7]);
    
%     %PLOT AXES LABELS
%     ha = handles.GUI.axes1;
%     h=gobjects(1,3);
%     Lr=input.cube.Lr;
%     h(1) = text(Lr(1),0,0,    ' FRONT','parent',ha); %1
%     %h(2) = text(-Lr(1),0,0,   ' BACK','parent',ha); %3
%     h(2) = text(0,Lr(2),0,    ' RIGHT','parent',ha); %2
%     %h(4) = text(0,-Lr(2),0,   ' LEFT','parent',ha); %4
%     h(3) = text(0,0,-Lr(3),   ' TOP','parent',ha); %5
%     %h(6) = text(0,0,Lr(3),    ' BOTTOM','parent',ha); %6
%     set(h,'FontSize',16,'Color',[1 1 1]*.6,'fontweight','normal');
%     handles.detectorfacelabels = h;
    
    %GET COLORS
    if G1.nutid>0
        c0=handles.GUI.figure1.Color;
        c=cell2mat(G1.upidcolor);
        pos=handles.GUI.promptColors.Position;  Cdata=ones(pos(4),pos(3),3);
        for i=1:3
            Cdata(:,:,i)=c0(i);
            for j=1:G1.nupid
                height = 17; %pixels
                pixelrows = 1+height*(j-1):height*j-1;
                Cdata(pixelrows,:,i) = c(j,i);
            end
        end
        handles.GUI.promptColors.CData = Cdata;
    end
    
    %PLOT UNIT VECTORS
    deleteh(handles.uvec);  handles.uvec=gobjects(1,G1.nua);
    p0 = G1.p0;    p = [p0; p0-G1.uvec(1,:)*70];
    handles.uvec(1) = plot3(p(:,1),p(:,2),p(:,3),'-','LineWidth',3,'Color',[.3 .3 .3]);  %original particle
    if G1.nua>0
        [~, c] = fcnpid2name(fcntid2pid(G1.ua,G1.tid,G1.pid));
        for i = 1:G1.nua
            p = [p0; p0+G1.uvec(i+1,:)]; %tid==1
            handles.uvec(i+1) = plot3(p(:,1),p(:,2),p(:,3),'-','LineWidth',1.5,'Color',c{i});
        end
    end
    
    %SET GUI CLIMS
    fcnsetGUIclims(input,output,handles,flags);
    h=legend(handles.GUI.axes1,handles.upid(isgraphics(handles.upid)),'Location','NorthEast');  %fcnlinewidth(h,5);  fcnmarkersize(h,16)
    h.EdgeColor=[1 1 1]*.8;  h.FontSize=12;
    h.Position(1:2) = [.90 .65];
    %handles.GUI.uitable1.Units='pixels';
    %handles.GUI.uitable1.Position(4) = (G1.nupid+1)*17 + 25;
    
    if ~isempty(findobj('type','figure','Name','Event Viewer'));  eventViewer(input,output); end
end

%FIND CURRENT PLOT TIME ---------------------------------------------------
input.plotTime = input.tMin + handles.GUI.slider1.Value*(input.tMax.current-input.tMin);
handles.GUI.text7.String = sprintf('t=%.4fns',input.plotTime);
handles.GUI.text4.String = sprintf('%.2fns',input.tMax.current);

%RUN THROUGH ALL THE BUTTONS ----------------------------------------------
if flags.update.pushbutton1==1;  flags.update.pushbutton1=0;  end
if flags.update.pushbutton2==1;  flags.update.pushbutton2=0;  end

%UPDATE PLOTS -------------------------------------------------------------
set(handles.upid(~G1.upidflag),'XData',[],'YData',[],'ZData',[]);
set(handles.upidsp,'XData',[],'YData',[],'ZData',[]);
set(handles.upidcp,'XData',[],'YData',[],'ZData',[]);

ptable = cell(G1.nupid+1,6);
photons.started = photons.sourceTime<input.plotTime;
for i = find(G1.upidflag)
    fcnUpdate1ParticleType(input, handle(handles.upid(i)), G1, i);
    ptable(i,:) = fcnPhotons(input,handles,flags,photons,PE, G1, i);
end
fcnupdateGUItable(G1,handles,ptable)

%DETECTOR OUTLINE
s = 'off';  if ischecked(handles.GUI.CCDMenudetectoroutline);  s='on';  end; handles.detectoroutline.Visible=s;

%PMTs
set(handles.gridanode,'XData',[],'YData',[],'ZData',[],'CData',[],'FaceVertexAlphaData',[]);
set(handles.stripanode,'XData',[],'YData',[],'ZData',[],'CData',[]);
if any(flags.status.plotside) && PE.n>0
    v1 = find( PE.t>(input.plotTime-.3-6*flags.status.lappd) & PE.t<input.plotTime);
    fcnplotpixels(input,output,handles,G1,flags,PE,v1);
end

%RESET SLIDER AND NEW EVENT FLAGS -----------------------------------------
flags.update.edit1=0;
end