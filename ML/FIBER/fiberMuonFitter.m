% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = fiberMuonFitter(input,output,flags,PE,G1,handles,x,s,plotflag)
results=[];
[~, si] = sortrows(x,4);  x=x(si,:);  s=s(si,:);  %sort by time
w=1./s(:,3).^2;

%VECTOR
[p0,vhat] = fitLineND(x(:,1:3),w);


%DIRECTION
p1 = fcndetectorintercept(input,p0,-vhat);
p2 = fcndetectorintercept(input,p0,vhat);
vhat = muonDirection(input,p0,vhat,p1,p2,x,s);


%ENERGY
rhat=norm(p2-p1); %mm
dEhat=sum(x(:,5)); %MeV


%ERROR
i = find(G1.pid==13,1,'first');
if any(i)
    v=fcnvec2uvec(G1.p2(i,:)-G1.p1(i,:));
    dE=sum(G1.de(G1.p1inside & G1.p2inside));
    
    p1=G1.p1(i,:);  vec = G1.p2(i,:)-p1;
    pn = fcndetectorintercept(input,p1,vec);  r = rangec(p1,pn);
    results.true = [v, dE, r];
    results.xhat = [vhat dEhat rhat];
end


%PLOTTING
if plotflag
    ct = vhat*v';  theta=acosd(ct);
    fcnplot3([p0-vhat*120; p0+vhat*120],'g-','linewidth',2.5);
    fprintf('%.1fdeg Muon error (%.1f el, %.1f az)\n',theta,fcnel(vhat)*r2d-fcnel(v)*r2d,fcnaz(vhat)*r2d-fcnaz(v)*r2d)
end

end



function vhat = muonDirection(input,p0,vhat,p1,p2,x,s)
% c = 299.792458; % (mm/ns)
% t=x(:,4);
% w=1./s(:,4).^2; w=w/sum(w);
% 
% that1 = rangec(p1,x(:,1:3))/c + min(t); %(ns)
% that2 = rangec(p2,x(:,1:3))/c + min(t); %(ns)
% 
% l1=weightedMean((that1-t).^2,w);
% l2=weightedMean((that2-t).^2,w);
% 
% if l1>l2; vhat=-vhat; end


n=size(x,1);  d=rangec(p1,x(:,1:3));  tv=(1:n)'.^2;

r1=weightedMean(d,tv);
r2=weightedMean(d,flipud(tv));
if r2>r1; vhat=-vhat; end
end