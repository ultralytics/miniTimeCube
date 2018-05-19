function A = oneBookend(F,I,D,PEs,t0,gain)
MTC766flag=false;
if MTC766flag
    [~,i] = fcnpruninglist;  D(:,~i)=0;  I(~i,2:end)=0;
end

A.D = D;
A.ZF = F;
A.ZI = I;

A.pid = F(:,1);
A.amplitude = F(:,2);
A.t = F(:,4);

A.A     = accumarray(F(:,1), F(:,2), [size(D,2) 1]); %Amplitude
%A.A     = max( accumarray(F(:,1), F(:,5), [size(D,2) 1]), 0); %Integral
A.N     = A.A./gain;
A.Ntrue = []; 
A.Nsum = numel(A.t);
A.Nsumtrue = sum(PEs);
A.t0 = t0;
