function [mu,s,k,X] = fcngmdistribution(X)
%fits data and tells you how many normals make up the underlying mixture distribution
mu=[];
s=[];
k=[];



n = 2; %max number points to try
BIC = inf(n,1);
obj = cell(n,1);
%fig; x = linspace(-200,200,100)'; fcnhist(X,x);
for k = 1:n
    f = gmdistribution.fit(X,k,'Replicates',5);
    %[x{nb}, fx] = patternsearch(@(x) fcnfermatpointvectorized(k,zt,x), x0, A, B ,Aeq, Beq, LB, UB, input.optimizer.psoptions2);
    f = gmdistribution.fit(rejectoutliers(f,X),k,'Replicates',5);
    
    BIC(k) = f.BIC;
    
    %fig; plot(x,pdf(f,x))
    obj{k}=f;
end
[~,k]=min(BIC)
f=obj{k};
mu = f.mu
s = f.Sigma;
X=rejectoutliers(f,X);
end



function X=rejectoutliers(f,X)
n = numel(f.PComponents);
idx = cluster(f,X);

for i=1:n
    j = find(idx==i);
    [~,inliers]=fcnsigmarejection(X(j),3,3);
    X(j(~inliers))=nan;
end
X=X(isfinite(X));
    
end