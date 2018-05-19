%INITIALIZATION FUNCTIONS -------------------------------------------------
function varargout = nView(varargin)
% nView M-file for nView.fig
%      nView, by itself, creates a new nView or raises the existing
%      singleton*.
%
%      H = nView returns the handle to a new nView or the handle to
%      the existing singleton*.
%
%      nView('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in nView.M with the given input arguments.
%
%      nView('Property','Value',...) creates a new nView or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nView_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nView_OpeningFcn via varargin.
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nView_OpeningFcn, ...
    'gui_OutputFcn',  @nView_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function nView_OpeningFcn(hObject, eventdata, handles, varargin)
fprintf('GEANT GUI Opening, Please Wait...\n');
disableButtons(handles)
handles.output = hObject;
guidata(hObject, handles);
function varargout = nView_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%set(handles.figure1,'Position',[10 50 1400 890])
%fprintf('Running nView_OutputFcn...\n')

%try
    pushbutton1_Callback([],[],handles)
% catch
%     fprintf('Failed to Open GUI (nview.m line 46) ... ')
%     close(gcf)
%     return
% end

%fprintf('GEANT GUI Opened.\n')

function figure1_CloseRequestFcn(hObject, eventdata, handles)
evalin('base','flags.status.play=0;') %stop rotation
evalin('base','delete(get(0,''Children'')); clear')
fprintf('GEANT GUI Closed.\n')


%FILE MENU FUNCTIONS ------------------------------------------------------
function OpenGEANTMenuItem_Callback(hObject, eventdata, h)
evalin('base','clear GEANT')
h.edit1.String='0'; %set eventnumber to 0
pushbutton1_Callback(hObject, eventdata, h) %reset

function OpenMenuItem_Callback(hObject, eventdata, h) %LOAD DETECTOR MODEL
path = evalin('base','input.directory');
[file, path] = uigetfile([path '/SAVED/*.mat']);
disableButtons(h)
if isequal(file, 0)
    enableButtons(h)
    return
else
    data = load([path file]);
end
input = data.input;
flags = data.flags;
input.eventNumber=0;

%FLUID
h.editFluid_yield.String = input.yield;
h.editFluid_ExponentialLambda.String = input.scintillatorDecay;
popupmenuFluid_Callback(hObject, false, h) %DO NOT RUN UPDATE VARS BECAUSE popupmenuPMT_Callback below runs it!

%CUBE
h.edit1.String = input.eventNumber;
h.edit51.String = sprintf('%.3g  %.3g  %.3g',input.cube.Lr*1E-3*2);
h.edit52.String = input.cube.requestedpixels;
h.editSensorModel_CoverageFraction.String=input.cube.coverageFraction;
h.popupmenuPMT.Value=find(strcmp(h.popupmenuPMT.String,input.cube.MCPname));
if strcmp('cylinder',input.cube.shape); h.ShapeCylinder.Checked='on'; else  h.ShapeCylinder.Checked='off'; end
if strcmp('cube',input.cube.shape);     h.ShapeCube.Checked='on'; else      h.ShapeCube.Checked='off'; end

if flags.status.plotside(1);     h.CCDMenu1.Checked='on'; else      h.CCDMenu1.Checked='off'; end
if flags.status.plotside(2);     h.CCDMenu2.Checked='on'; else      h.CCDMenu2.Checked='off'; end
if flags.status.plotside(3);     h.CCDMenu3.Checked='on'; else      h.CCDMenu3.Checked='off'; end
if flags.status.plotside(4);     h.CCDMenu4.Checked='on'; else      h.CCDMenu4.Checked='off'; end
if flags.status.plotside(5);     h.CCDMenu5.Checked='on'; else      h.CCDMenu5.Checked='off'; end
if flags.status.plotside(6);     h.CCDMenu6.Checked='on'; else      h.CCDMenu6.Checked='off'; end
if flags.status.plotside(6);     h.CCDMenu6.Checked='on'; else      h.CCDMenu6.Checked='off'; end
if isfield(flags.status,'fibers') && flags.status.fibers;          h.FiberFlag.Checked='on'; else     h.FiberFlag.Checked='off'; end

popupmenuPMT_Callback(hObject, eventdata, h)
enableButtons(h)
function SaveMenuItem_Callback(hObject, eventdata, h)
disableButtons(h)
evalin('base','fcnsavefile(input,flags,handles);')
enableButtons(h)
function CCDMenu1_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenu2_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenu3_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenu4_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenu5_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenu6_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenudetectoroutline_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function CCDMenuAll_Callback(hObject, eventdata, handles)
if ischecked(hObject)
    set(hObject,'Checked','off')
    set(handles.CCDMenu1,'Checked','off')
    set(handles.CCDMenu2,'Checked','off')
    set(handles.CCDMenu3,'Checked','off')
    set(handles.CCDMenu4,'Checked','off')
    set(handles.CCDMenu5,'Checked','off')
    set(handles.CCDMenu6,'Checked','off')
else
    set(hObject,'Checked','on')
    set(handles.CCDMenu1,'Checked','on')
    set(handles.CCDMenu2,'Checked','on')
    set(handles.CCDMenu3,'Checked','on')
    set(handles.CCDMenu4,'Checked','on')
    set(handles.CCDMenu5,'Checked','on')
    set(handles.CCDMenu6,'Checked','on')
end
updateVars(handles)
function CCDMenuSolidAngle_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function PhotonsMenuSize4_Callback(hObject, eventdata, handles)
set(handles.PhotonsMenuSize4,'Checked','off')
set(handles.PhotonsMenuSize8,'Checked','off')
checkuncheckUpdate(hObject,handles)
function PhotonsMenuSize8_Callback(hObject, eventdata, handles)
set(handles.PhotonsMenuSize4,'Checked','off')
set(handles.PhotonsMenuSize8,'Checked','off')
checkuncheckUpdate(hObject,handles)
function PhotonMenuShowAlive_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function PhotonMenuShowSource_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function PhotonMenuShowEnd_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
function PhotonMenuColorBy_Callback(hObject, eventdata, h)
set(h.PhotonMenuColorBySource,'Checked','off')
set(h.PhotonMenuColorByGen,'Checked','off')
set(h.photonMenuColorByWavelength,'Checked','off')
set(h.PhotonMenuColorByAge,'Checked','off')
checkuncheckUpdate(hObject,h)
function checkuncheckUpdate(hObject,handles)
fcncheckuncheck(hObject);
updateVars(handles)
function AxisOn_Callback(hObject, eventdata, handles)
if ischecked(hObject)
    axis off
else
    axis on
end
fcncheckuncheck(hObject)


function realdataflag_Callback(hObject, eventdata, handles)
evalin('base','flags.update.edit1=1;')
checkuncheckUpdate(hObject,handles)
function realtimevoltage_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)


%CREATE GUI FUNCTIONS -----------------------------------------------------
function figure1_CreateFcn(hObject, eventdata, h)
function popupmenu1_CreateFcn(hObject, eventdata, handles)
function slider1_CreateFcn(hObject, eventdata, handles)
hObject.Value=1;
function edit1_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.eventNumber=' hObject.String ';'])
function edit51_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.cube.Lr=[' hObject.String ']*1E3/2;'])
function edit52_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.cube.requestedpixels=' hObject.String ';'])
function editFluid_ExponentialLambda_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.scintillatorDecay=' hObject.String ';'])
function editFluid_yield_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.yield=' hObject.String ';'])
function editSensorModel_CoverageFraction_CreateFcn(hObject, eventdata, handles)
evalin('base',['input.cube.coverageFraction=' hObject.String ';'])



%CALLBACK FUNCTIONS -------------------------------------------------------
function pushbutton1_Callback(hObject, eventdata, h) %RESET
rng(1); cla; clc

%PASS ALL HANDLES TO WORKSPACE
assignin('base','hgui',h); evalin('base','handles.GUI=hgui; clear hgui;')
evalin('base','nViewInit')
updateVars(h)
setDt(h)
function pushbutton3_Callback(hObject, eventdata, handles) %AXES TIGHT
disableButtons(handles)
aLims1=axis; %current axis;
aLims1(5:6) = get(handles.axes1,'ZLim');
if strcmp(hObject.String,'Axis Tight')% zoom out
    axis(handles.axes1,'image')
    aLims2 = axis; %desired axis
    aLims2(5:6) = get(handles.axes1,'ZLim');
    set(hObject,'String','Axis Cube')
    for n = logspace(0,-3,20) %interpolate between the two "zoom in" or "zoom out"
        axis(handles.axes1, aLims1*(n)+aLims2*(1-n))
        pause(.01)
    end
else %zoom in
    cp = evalin('base','input.cube.Lr;');
    aLims2=[-1 1 -1 1 -1 1]*10 + [cp(1)*[-1 1] cp(2)*[-1 1] cp(3)*[-1 1]]; %desired axis
    set(hObject,'String','Axis Tight')
    for n = logspace(0,-3,20) %interpolate between the two "zoom in" or "zoom out"
        axis(handles.axes1, aLims1*(n)+aLims2*(1-n))
        pause(.01)
    end
end
enableButtons(handles)
function pushbuttonRotateLeft_Callback(~, ~, ~) %#ok<*DEFNU>
[az,el] = view;  fcnrotateaxes(gca,[az+90 el], 30)
function pushbuttonRotateRight_Callback(~, ~, ~)
[az,el] = view;  fcnrotateaxes(gca,[az-90 el], 30)
function pushbuttonRotateTop_Callback(~, ~, ~)
[az,el] = view;  fcnrotateaxes(gca,[az min(el+45,90)], 25)
function pushbuttonRotateBottom_Callback(~, ~, ~)
[az,el] = view;  fcnrotateaxes(gca,[az max(el-45,-90)], 25)
function popupmenu1_Callback(hObject, eventdata, handles) %SELECT MEV
disableButtons(handles)
evalin('base',['GEANT=fcnloadGEANT(GEANT.files.names{' num2str(get(handles.popupmenu1, 'Value')) '});'])
edit1_Callback(hObject, eventdata, handles)
function slider1_Callback(hObject, eventdata, handles)
evalin('base','[input,output,G1,photons,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles); fcnResizeAxis(flags,handles,input,G1)')
%updateVars(handles)
function edit1_Callback(hObject, eventdata, handles)
evalin('base','flags.update.edit1=1;')
updateVars(handles)
function edit51_Callback(hObject, eventdata, handles)
popupmenuPMT_Callback(hObject, eventdata, handles)
function edit52_Callback(hObject, eventdata, handles)
popupmenuPMT_Callback(hObject, eventdata, handles)
function editFluid_ExponentialLambda_Callback(hObject, eventdata, handles)
popupmenuFluid_Callback(hObject, eventdata, handles)
evalin('base','fcnupdateAxes2(input,handles,G1,photons);')
function editFluid_yield_Callback(hObject, eventdata, handles)
popupmenuFluid_Callback(hObject, eventdata, handles)
function editSensorModel_CoverageFraction_Callback(hObject, eventdata, handles)
popupmenuPMT_Callback(hObject, eventdata, handles)

%CHECKBOX CALLBACK FUNCTIONS ----------------------------------------------
function checkbox52_Callback(hObject, eventdata, handles) %VERTEX ON WALL
edit1_Callback([],[],handles)
function checkbox58_Callback(hObject, eventdata, handles) %VERTEX CENTERED
edit1_Callback([],[],handles)
function checkbox59_Callback(hObject, eventdata, handles) %VERTEX FIXED ANGLE
edit1_Callback([],[],handles)

function checkboxautoresize_Callback(hObject, eventdata, handles) %AUTO RESIZE
disableButtons(handles)
set(handles.axes1,'CameraPositionMode','auto')
set(handles.axes1,'CameraTargetMode','auto')
set(handles.pushbutton3,'String','Axis Even')
evalin('base','fcnResizeAxis(flags,handles,input,G1)')
enableButtons(handles)
function checkboxPlay_Callback(hObject, eventdata, handles)
a = load('GUIplaypauseCdata.mat');
if get(hObject,'Value')
    set(hObject,'Cdata',a.pause);
    evalin('base','flags.status.play=1; play')
    set(hObject,'Cdata',a.play);
else
    evalin('base','flags.status.play=0;')
    set(hObject,'Cdata',a.play);
end
function checkboxDtDown_Callback(ho, eventdata, h)
setDt(h,h.textDt.Value/2)
function checkboxDtUp_Callback(ho, eventdata, h)
setDt(h,h.textDt.Value*2)
function setDt(h,dt)
if nargin==1 %set workspace dt to the current value if none passed in
    dt = h.textDt.Value;
end
set(h.textDt,'Value',dt,'String',['dt=' num2str(dt) 'ns']);
evalin('base',['input.dt=' num2str(dt) ';'])
%LOCAL FUNCTION ONLY ------------------------------------------------------
function disableButtons(handles)
% handles.axes1.HitTest='off';
% handles.figure1.Pointer='watch';
% handles.pushbutton1.Enable='off'; %reset?
% handles.pushbutton3.Enable='off'; 
% handles.slider1.Enable='off';
% handles.popupmenu1.Enable='off';
% handles.edit1.Enable='off';
% handles.edit51.Enable='off';
% handles.edit52.Enable='off';
% handles.editFluid_ExponentialLambda.Enable='off';
% handles.editFluid_yield.Enable='off';
% handles.editSensorModel_CoverageFraction.Enable='off';
% handles.checkboxautoresize.Enable='off';
%set(findobj(gcf,'Type','uicontrol'),'Enable','off')
%drawnow
function enableButtons(handles)
% handles.axes1.HitTest='on';
% handles.figure1.Pointer='arrow';
% handles.pushbutton1.Enable='on'; %reset?
% handles.pushbutton3.Enable='on'; 
% handles.slider1.Enable='on';
% handles.popupmenu1.Enable='on';
% handles.edit1.Enable='on';
% handles.edit51.Enable='on';
% handles.edit52.Enable='on';
% handles.editFluid_ExponentialLambda.Enable='on';
% handles.editFluid_yield.Enable='on';
% handles.editSensorModel_CoverageFraction.Enable='on';
% handles.checkboxautoresize.Enable='on';
%findobj(gcf,'Type','uicontrol'),'Enable','on')
function updateVars(handles)
disableButtons(handles)

mmen = [0 1E6];
if ischecked(handles.realdataflag)==0
    mmen = evalin('base','minmax3(GEANT.ve)'); %minmax eventnumber
    handles.text8.String=sprintf('Event (%g-%g)',mmen);
end
str = num2str( floorandceil(eval(handles.edit1.String),mmen(1),mmen(2)) ); %max min
handles.edit1.String=str;
evalin('base',['input.eventNumber=' str ';'])


evalin('base',['flags.status.vertexonwall=' num2str(handles.checkbox52.Value) ';']) %vertex on wall
evalin('base',['flags.status.vertexcentered=' num2str(handles.checkbox58.Value) ';']) %vertex centered
evalin('base',['flags.status.vertexfixedangle=' num2str(handles.checkbox59.Value) ';']) %vertex fixed angle

evalin('base',['flags.status.plotcurrentvoltage=' num2str(ischecked(handles.realtimevoltage)) ';'])

evalin('base',['flags.status.plotside(1)=' num2str(ischecked(handles.CCDMenu1)) ';'])
evalin('base',['flags.status.plotside(2)=' num2str(ischecked(handles.CCDMenu2)) ';'])
evalin('base',['flags.status.plotside(3)=' num2str(ischecked(handles.CCDMenu3)) ';'])
evalin('base',['flags.status.plotside(4)=' num2str(ischecked(handles.CCDMenu4)) ';'])
evalin('base',['flags.status.plotside(5)=' num2str(ischecked(handles.CCDMenu5)) ';'])
evalin('base',['flags.status.plotside(6)=' num2str(ischecked(handles.CCDMenu6)) ';'])
evalin('base',['flags.status.plotdetectoroutline=' num2str(ischecked(handles.CCDMenudetectoroutline)) ';'])
evalin('base',['flags.status.showSolidAngles=' num2str(ischecked(handles.CCDMenuSolidAngle)) ';'])

evalin('base',['flags.status.plotphotonsalive=' num2str(ischecked(handles.PhotonMenuShowAlive)) ';'])
evalin('base',['flags.status.plotphotonsource=' num2str(ischecked(handles.PhotonMenuShowSource)) ';'])
evalin('base',['flags.status.plotphotonend=' num2str(ischecked(handles.PhotonMenuShowEnd)) ';'])
evalin('base',['flags.status.photoncolorbysource=' num2str(ischecked(handles.PhotonMenuColorBySource)) ';'])
evalin('base',['flags.status.photoncolorbygen=' num2str(ischecked(handles.PhotonMenuColorByGen)) ';'])
evalin('base',['flags.status.photoncolorbywavelength=' num2str(ischecked(handles.photonMenuColorByWavelength)) ';'])
evalin('base',['flags.status.photoncolorbyage=' num2str(ischecked(handles.PhotonMenuColorByAge)) ';'])

evalin('base',['if ' num2str(ischecked(handles.PhotonsMenuSize4)) '; input.photonsize=4; end'])
evalin('base',['if ' num2str(ischecked(handles.PhotonsMenuSize8)) '; input.photonsize=8; end'])

cubeLength = handles.edit51.String;  clv = eval(['[' cubeLength ']']); %cube length vector
if numel(clv)~=3; %someone input cube dimensions of '1' instead of '1 1 1'
   clv = mean(clv)*[1 1 1];  
   cubeLength = sprintf('%.3g  %.3g  %.3g',clv);
   handles.edit51.String=cubeLength;
end
evalin('base',['input.cube.Lr=[' cubeLength ']*1E3/2;'])

evalin('base',['input.cube.requestedpixels=' handles.edit52.String ';'])
evalin('base',['input.cube.coverageFraction=' handles.editSensorModel_CoverageFraction.String ';'])
evalin('base',['input.yield=' handles.editFluid_yield.String ';'])
evalin('base',['input.scintillatorDecay=' handles.editFluid_ExponentialLambda.String ';'])

evalin('base',sprintf('input.volumeFilename=''%s'';', strtrim(handles.popupmenuFluid.String{handles.popupmenuFluid.Value})) )
evalin('base','input=fcninittables(input); flags.update.edit1=1;')

evalin('base','fcnsetGUIclims(input,output,handles,flags); [input,output,G1,photons,PE,flags,handles]=fcnupdateAxes1(input,output,GEANT,G1,photons,PE,flags,handles); fcnupdateAxes2(input,handles,G1,photons); fcnResizeAxis(flags,handles,input,G1)')
handles.edit52.String = evalin('base','input.cube.pixels'); %adjust pixel count
enableButtons(handles)



%UNDER CONSTRUCTION -------------------------------------------------------
function figure1_WindowScrollWheelFcn(hobj, event, h) 
% factor = event.VerticalScrollCount;
% x = get(h.axes1,'CurrentPoint');
% set(h.axes1,'CameraTargetMode','manual','CameraTarget',mean(x,1))
% if factor>0
%     amount = 1.4;
% else
%     amount= 0.6;
% end
% zoom(amount)
% centerPointer(h)
% function centerPointer(h)
% set(h.figure1,'Units','normalized'); p2=get(h.figure1,'Position');
% set(h.axes1,'Units','normalized'); p1 = get(h.axes1,'Position');
% ss=get(0,'ScreenSize');
% p2 = p2.*ss([3 4 3 4]);
% p1 = p1.*p2([3 4 3 4]);
% p3 = [p1(1)+p1(3)/2, p1(2)+p1(4)/2] + p2(1:2);
% set(0,'PointerLocation',p3);


function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
'moving'


%ADDCHECKBOX CHECKLIST
%1. add to enableButtons and disableButtons callbacks
%2. add to updateVars callback
%4. add to clearInvalidCheckboxes callback
%6. add to pushbutton1
%7. add to init.m


% --------------------------------------------------------------------
% --------------------------------------------------------------------


function popupmenuFluid_Callback(hObject, eventdata, handles)
updateVars(handles);

function popupmenuFluid_CreateFcn(hObject, eventdata, handles)
s = what([fcnpathm1(mfilename('fullpath')) 'LIBRARIES/MATERIALS/']); rownames=s.mat; rownames{end+1}='manual';
hObject.String=rownames;
hObject.Value=1; %pick first row
evalin('base',sprintf('input.volumeFilename=''%s'';', strtrim(rownames{hObject.Value})) )

function popupmenuPMT_Callback(hObject, eventdata, handles)
cells = handles.popupmenuPMT.String;
evalin('base',sprintf('input.cube.MCPname=''%s'';', strtrim(cells{handles.popupmenuPMT.Value})) )
evalin('base','input= fcninittables(input); flags.update.edit1=1; flags.update.detectorgeometry=1;')
updateVars(handles);  

function popupmenuPMT_CreateFcn(hObject, eventdata, handles)
s = what([fcnpathm1(mfilename('fullpath')) 'LIBRARIES/PMTS/']); rownames=s.mat; rownames{end+1}='manual';
hObject.String=rownames;
hObject.Value=numel(rownames); %pick last row
evalin('base',sprintf('input.cube.MCPname=''%s'';', strtrim(rownames{hObject.Value})) )

function promptColors_Callback(hObject, eventdata, handles)


%NEW TABLE!!!!!!!!
% --------------------------------------------------------------------
function uitable1_CellSelectionCallback(h, eventdata, handles)
[nr,nc] = size(get(h,'Data'));
if ~isempty(eventdata.Indices)
    row=eventdata.Indices(1);
    if row<nr
        row0=find(evalin('base',sprintf('G1.upid==G1.upid(%.0f)',row)));
        evalin('base',sprintf('G1.upidflag(%.0f)=~G1.upidflag(%.0f);',row0,row0))
        updateVars(handles)
    end
end


function pushbuttonnextevent_Callback(hObject, eventdata, handles)
n = eval(handles.edit1.String);
handles.edit1.String=num2str(n+1);
edit1_Callback(hObject, eventdata, handles)


function ShapeCube_Callback(hObject, ~, h)
h.ShapeCube.Checked = 'on';
h.ShapeCylinder.Checked = 'off';
popupmenuPMT_Callback([], [], h)


function ShapeCylinder_Callback(hObject, ~, h)
h.ShapeCube.Checked = 'off';
h.ShapeCylinder.Checked = 'on';
popupmenuPMT_Callback([], [], h)


% --------------------------------------------------------------------
function FiberFlag_Callback(hObject, eventdata, handles)
checkuncheckUpdate(hObject,handles)
