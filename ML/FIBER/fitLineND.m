function [mu,V] = fitLineND(x,w)
%x = [n x nd] vector, nd=number dimensions
if numel(x(:,1))==1
    fprintf('WARNING: fitLineND.m being asked to fit a line to 1 point\n')
    mu = x;
    V = [0 0 0];
    return
end

C = weightedcov(x,w);
[eigvec,eigval] = eigs(C);

mu=weightedMean(x,w,1);
%plot3(mu(1),mu(2),mu(3),'g.','markersize',50);

[~,i]=max(diag(eigval));
V = eigvec(:,i)';

% %PLOT
% A=mu-V*120;
% B=mu+V*120;
% fcnplot3([A; B],'g-','linewidth',2.5);


