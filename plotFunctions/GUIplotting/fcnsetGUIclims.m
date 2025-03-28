% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnsetGUIclims(input,output,handles,flags)
if isempty(output) || ~isfield(output,'D'); return; end

MTCflag = ischecked(handles.GUI.realdataflag);

if MTCflag
    handles.GUI.axes1.CLim = minmax3(output(1).D(:)) + [0 eps]; %volts
else
    handles.GUI.axes1.CLim = [0 input.cube.dsp.dynamicrange(2)]; %volts
end

