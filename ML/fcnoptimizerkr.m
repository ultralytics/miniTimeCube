% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function k = fcnoptimizerkr(input,output)
%reduced version of optimizer constants using only active pixels
k = fcnoptimizerk(input,output(1).N);

ti = output.t~=0;  
k.zt = output.t(ti);
[i, ~, j] = fcnunique(output.pid(ti));


k.upid = i; %unique pixel id's
k.nhpp = output.N(i)'; %number of hits per unique pixel
k.pid = uint32(j); %pixel id's

k.pxyz = k.pxyz(i,:);
k.nxyz = k.nxyz(i,:);

k.npixels = numel(i);
k.nZ = numel(j);
k.ovnz = ones(k.nZ,1,'uint8');

k.QEmap = k.QEmap(i);
