% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [x, fx] = fcnMLenergy(input,p0,zN,x0) %vectorized by default
if sum(zN>0)<3; x=0; fx=0; return; end
k = fcnoptimizerk(input,zN);  k.reflections = 0;

n = size(p0,1);
if nargin==3
    x0 = zeros(1,n);
    for i=1:n
        x0(i) = fcnguessinitenergy(k,zN,p0(i,:));
    end
end

LB=.001*ones(1,n);  UB=15*ones(1,n);
[x, fx] = patternsearch(@(x) fcnL(k,p0,zN,x), x0, [], [], [], [], LB, UB, [], input.optimizer.psoptions3);
end


function fx = fcnL(k,p,zN,X)
%X = [xyze, xyze ...];
[nv,np] = size(X);  np=np/1; %nv = number of vectorized points
iv      = (0:np-1)'*1;

[~, I, J] = fcnunique(p(:,1:3),'rows');

F0 = fcnsolidanglevectorized(k,p(I,:));
F0 = sum(F0,3); %sum any reflections

e       = fcnrow(X(:,iv+1))*k.yield;
wf      = F0(:,J).*e;
lambda  = sum(reshape(wf,[numel(zN) nv np]),3) + 1E-323;

fx = -sum( zN.*log(lambda) - (lambda+k.poiss.gammaln) ); %NLL from MATLAB poisspdf function. SLOWER, WORKS WITH NON INTEGER ZN's
%fx = -( sum(log(lambda(k.poiss.z1,:))) + k.poiss.zN2'*log(lambda(k.poiss.z2,:)) - sum(lambda) - k.poiss.gammalnsum ); %NLL FASTER
end


%EXTENDED POISSON ---------------------------------------------------------
% fig; lambda = 2.5;
% x=0:1:6;  stem(x,poisspdf(x,lambda),'b.','markersize',30)
% zN=0:.01:6;  y=exp( (zN.*log(lambda) - lambda - gammaln(zN+1))); plot(zN,y,'g-'); 
% fcnlinewidth(2);  legend('Poisson Distribution,  f(k,  \lambda) = \lambda^ke^{-\lambda} / k!','Extended Poisson Distribution,  f(k,  \lambda) = e^{k log(\lambda) - \lambda - log\Gamma(k+1)}')
% xyzlabel('k','pdf','',sprintf('Extended Poisson Example  \\lambda = %.1f',lambda)); fcnfontsize(16);