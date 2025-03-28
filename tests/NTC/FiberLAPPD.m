% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = FiberLAPPD(input,output,handles,PE,flags,plotflag)
if plotflag; closeallexcept(handles.GUI.figure1);  deleteh(findobj(handles.GUI.figure1,'DisplayName','test')); end
A=output(1).D;
[nx, ns, ~, nlappd] = size(A); %number of times, number of strips (30), number sides (2), number of lappd's

%PULSE RECONSTRUCTION
%B=fcnSplinePulseReconstruction(reshape(A,[nx ns*2*nlappd]),input.cube.dsp.dynamicrange(2));  A=reshape(B,[nx ns 2 nlappd]);

end

