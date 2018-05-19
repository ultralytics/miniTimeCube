function [] = fcnupdateGUItable(G1,handles,x)
h = handles.GUI.uitable1;  if G1.nupid==0; h.Data=[]; return; end

s = sum(cell2mat(x),1); if isempty(s); s=zeros(1,5); end
x(1:end-1,1) = G1.upidname; %names
x(end,:) = {'Sum',s(1),s(2),s(3),s(4),s(5)};

h.Data=x;
%set(h,'ColumnName',{'Name','TID','PID','dE(MeV)','Photons','PEs'})
%set(handles.GUI.editFluid_RefractiveIndex,'Selected','off')
