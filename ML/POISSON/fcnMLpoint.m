function [x, fx] = fcnMLpoint(input,output,x0) %vectorized by default
k = fcnoptimizerk(input,output);  k.reflections = 1;
zN = k.nhpp;

if nargin==2
    s = sqrt(zN);
    p0 = s'*k.pxyz/sum(s) * .95 * 0;
    e0 = fcnguessinitenergy(k,zN,p0);
    x0 = [p0 e0];
end

LB=[]; UB=[];
[x, fx]  = patternsearch(@(x) fcnL(k,zN,x), x0, [], [], [], [], LB, UB, @(x) nonlcon(x,input), input.optimizer.psoptions3);
%[x, fx] =       fmincon(@(x) fcnL(k,zN,x), x0, [], [], [], [], LB, UB, @(x) nonlcon(x,input), input.optimizer.options2);
end


function fx = fcnL(k,zN,X)
%X = [xyze, xyze ...];
[nv,np] = size(X);  np=np/4; %nv = number of vectorized points
iv      = (0:np-1)'*4;
posi    = (1:3)+iv;

p = reshape(X(:,posi),[np*nv 3]);
[~, I, J] = fcnunique(p(:,1:3),'rows');
Fs = fcnsolidanglevectorized(k,p(I,:));
Fs = sum(Fs,3); %sum any reflections

e       = fcnrow(X(:,iv+4))*k.yield;
wf      = Fs(:,J).*e;
lambda  = sum(reshape(wf,[numel(zN) nv np]),3) + 1E-323;

fx = -sum( zN.*log(lambda) - (lambda+k.poiss.gammaln) ); %NLL from MATLAB poisspdf function. SLOWER, WORKS WITH NON INTEGER ZN's
end


function [c,ceq]=nonlcon(x,input)
es = 1; %(mm) airgap between xhat and detector wall
if input.cube.shapeID==3 %cylinder
    c(1) = norm(x(1:2)) - input.cube.radius+es;
    c(2) = abs(x(3)) - input.cube.height/2+es;
    c(3) = .001 - x(4); %1keV floor
    c(4) = -50 + x(4); %50mev ceiling
else
    c = zeros(size(x,1),5);
    Lr = input.cube.Lr-es;
    c(:,1) = abs(x(:,1)) - Lr(1);
    c(:,2) = abs(x(:,2)) - Lr(2);
    c(:,3) = abs(x(:,3)) - Lr(3);
    c(:,4) = .001 - x(:,4); %1keV floor
    c(:,5) = -50 + x(:,4); %50mev ceiling
end
ceq = [];
end