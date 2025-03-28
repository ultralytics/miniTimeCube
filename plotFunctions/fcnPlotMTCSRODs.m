% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnPlotMTCSRODs(input,output,flags,photons,handles)
sa = {'SCROD','Row','Column','Channel','PMT','PMTRow','PMTCol','ASIC','ASICd','SOG'};
[SCROD, Row, Column, Channel, PMT, PMTRow, PMTCol, ASIC, ASICd, ~, SOG] = MTCpixelID2SRCCH(1:1536); %#ok<ASGLU>

load('MTCmapping.mat');  X=sortrows(X,10);  vSCROD=unique(X(:,2),'stable');  %vSCROD  = [1 2 4 6 7 8 10 11 12 255 101 0];
vSCROD      = [1 2 4 6 7 8 10 11 12 255 101 0];
vRow        = [0 1 2 3];
vColumn     = [0 1 2 3];
vChannel    = [0 1 2 3 4 5 6 7]; %#ok<*NASGU>
vPMT        = 1:24;
vPMTRow     = 1:8;
vPMTCol     = 1:8;
vASIC       = 1:16;
vASICd      = 1:192;
vSOG        = 1:3; %SCROD offset group

[ha,hf] = fig(3,4,1.6,1.2);
popoutsubplot(handles.GUI.axes1,ha(1))
%if isfield(input,'MTC'); fname=input.MTC.filename; else fname='GEANT'; end
fname = '';
title(sprintf('''%s'' event %f',fname,input.eventNumber))


sca(ha(2)) %color by unique channel
popoutsubplot(handles.GUI.axes1,gca); cla;
Cdata = 1:1536; %unique pixel number
hi = patch(input.cube.all.v.x', ...
    input.cube.all.v.y', ...
    input.cube.all.v.z','w', ...
    'EdgeColor',[.8 .8 .8], ... 'EdgeAlpha',0, ...
    'FaceColor','Flat', ...
    'FaceAlpha','Flat', ...
    'FaceVertexCData',Cdata', ...
    'FaceVertexAlphaData',.7, ...
    'CDataMapping','scaled', ...
    'AlphaDataMapping','none', ...
    'FaceLighting','none', ...
    'BackFaceLighting','unlit');  set(gca,'clim',[1 1536]);
hc=colorbar('East'); %fcncolorbar(.7,' pixel ',hc);
title('unique pixel IDs')


for pi=1:numel(sa)
    sca(ha(pi+2));
    %popoutsubplot(handles.GUI.axes1,gca); cla
    
    sb = sa{pi};
    eval(sprintf('a=%s; va=v%s;',sb,sb));  ua = unique(a);  nua=numel(ua);  nva=numel(va);
    for i = 1:nua
        j = find(a==ua(i));  k = find(va==ua(i));
        
        if any(j) && any(k)
            hi = patch(input.cube.all.v.x(j,:)', ...
                input.cube.all.v.y(j,:)', ...
                input.cube.all.v.z(j,:)','w', ...
                'EdgeColor',[.8 .8 .8], ... 'EdgeAlpha',0, ...
                'FaceColor',fcndefaultcolors(k,nva), ...
                'FaceAlpha','Flat', ...
                'FaceVertexAlphaData',.7, ...
                'AlphaDataMapping','none', ...
                'FaceLighting','none', ...
                'BackFaceLighting','unlit');
            set(hi,'DisplayName',sprintf('%s %g',sb,va(k)));
        end
    end
    axis equal 
    hli = fcnlegend(gca,1,'unique');
    set(hli,'Position',get(hli,'Position') - [.025 0 0 0])
    legend boxoff
end
fcnmatchlims(ha(1),ha(2:end))

set(ha,'cameraviewangle',7);
fcnfontsize(9)
end
