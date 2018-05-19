legend hide
% %START PLAYBACK WHEN NEUTRINO ENTERS DETECTOR
% input.plotTime = get(handles.GUI.slider1,'Value')*input.tMax.current;
% if input.plotTime == 0;
%     c = 299.792458; %mm/ns
%     input.tMin = G1.p0(1) - max(input.cube.Lr)/c;
%     input.plotTime = input.tMin;
% end

vidflag = false;
rotateflag = false;
fourflag = false; %all 4 viewing angles
if vidflag; clc
    hc=get(0,'children'); close(hc(hc~=handles.GUI.figure1));
    %fname=fcnincrementfname('Detector Model.mp4'); vidObj = VideoWriter(fname,'MPEG-4');  set(vidObj,'Quality',100);  open(vidObj);
    fname=fcnincrementfname('Detector Video.'); vidObj = VideoWriter(fname,'MPEG-4');  open(vidObj);
    import java.awt.Robot;  mouse=Robot;
    if fourflag
        [ha, hf] = fig(2,2,2.2,1.5); %colorbar('east');
    else
        [ha, hf] = fig(1,1,2.2,1.5); %colorbar('east');
    end
    set(hf,'Units','pixels');
    set(hf,'Position',[10 50 1920 1080]);
    axes(handles.GUI.axes1); firstframeflag = true;
    fi = 0;
end


if input.plotTime==input.tMax.current; input.plotTime=0; set(handles.GUI.slider1,'Value',0); end %restart if someone hits play and the marker is at the very end
tic
while flags.status.play% && input.plotTime<2
    input.plotTime = input.plotTime + input.dt;
    if input.plotTime<0
        set(handles.GUI.slider1,'Value',(input.plotTime-input.tMin)/(input.tMax.current-input.tMin))
        [input,~,~,~,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles);
    else
        if input.plotTime<(input.tMax.current + input.dt)
            sliderValue = min(get(handles.GUI.slider1,'Value'), 1);
            input.plotTime = input.tMin + sliderValue*(input.tMax.current-input.tMin); %get current spot if user moves it!
            input.plotTime = min( input.plotTime + input.dt, input.tMax.current);
            
            sliderValue = (input.plotTime-input.tMin)/(input.tMax.current-input.tMin);
            set(handles.GUI.slider1,'Value',sliderValue)
            [input,~,~,~,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles);
        else %reached end of slider!
            input.plotTime = input.tMax.current;
            set(handles.GUI.slider1,'Value',1)
            [input,~,~,~,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles);
            drawnow
            break
        end
    end

    %ROTATE PLOT
    if rotateflag %rotate plot
        [az,el] = view(handles.GUI.axes1);  view(handles.GUI.axes1,az-.5,el+.1);
    end
    %drawnow
    pause(.01)
    
%     if input.plotTime==input.dt || any( abs(input.plotTime-(0:.2:4)) < .001);
%         popoutsubplot(handles.GUI.axes1, fig)    
%         title(sprintf('t = %.3fns',input.plotTime)); colorbar
%         export_fig(gcf,'-q95','-r150','-a4',fcnincrementfname('exported.jpg'))
%         close(gcf)
%     end

    if vidflag
        %frame = getscreen(get(handles.GUI.figure1,'position') + [0 0 0 80]);  writeVideo(vidObj,frame);
        textc = min(fi/100,.7)*[1 1 1]; %fade text black to gray
        for i=1:numel(ha); cla(ha(i)); end
                
        popoutsubplot(handles.GUI.axes1, ha(1)); legend off
        if fourflag
            popoutsubplot(handles.GUI.axes1, ha(2)); view(ha(2),-90,90); legend off;    hts(1) = text(0,0,-2,'     TOP','color',textc,'fontsize',80,'parent',ha(2));
            popoutsubplot(handles.GUI.axes1, ha(3)); view(ha(3),-180,0); legend off;    hts(2) = text(0,0,-2,'     LEFT','color',textc,'fontsize',80,'parent',ha(3));
            popoutsubplot(handles.GUI.axes1, ha(4)); view(ha(4),-90,0); legend off;     hts(3) = text(0,0,-2,'     BACK','color',textc,'fontsize',80,'parent',ha(4));
        end
        text(0,0,-90,sprintf('%2.3fns         ',input.plotTime),'color',textc,'fontsize',80,'parent',ha(1));
        if firstframeflag; 
            set(ha,'cameraviewangle',get(ha(1),'cameraviewangle')*.65);
            hcb=colorbar('peer',ha(1),'west','FontSize',16); fcncolorbar(.7,'V',hcb)
            %colormap(hvid(1),hsv); 
            axis(ha,'off')
            firstframeflag=false; 
        end
        fcncolorbar(.7,[],hcb)
        %[az, el] = view(ha(1)); fcnrotateaxes(ha(1),[az+1 el],1);
        writeVideo(vidObj,getframe(hf));
        %figure(handles.GUI.figure1)
        %sca(handles.GUI.axes1)
        if input.plotTime>35; break; end
        fi = fi+1; if mod(fi,100)==0;  mouse.mouseMove(round(rand*100)+1,round(rand*100)+1); end %move mouse!
    end
    
    %if input.plotTime>5; break; end
end
if exist('vidObj','var');  close(vidObj); end %#ok<*UNRCH>
toc
legend show


input.tMin = 0;
input.plotTime = max(0,input.plotTime);
set(handles.GUI.slider1,'Value',(input.plotTime-input.tMin)/(input.tMax.current-input.tMin))
[input,output,G1,photons,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles);

fcnResizeAxis(flags,handles,input,G1)
 
clear az el tNow t vidObj vidflag sliderValue ha hf fi mouse fourflag rotateflag