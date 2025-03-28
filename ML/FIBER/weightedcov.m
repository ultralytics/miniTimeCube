% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function C = weightedcov(x,w)
%x = NxM, weights = Nx1
[N,M]=size(x);

%x = x-mean(x);;           C = x'*x/(N-1); %unweighted cov
x = x-weightedMean(x,w);   C = sum(w)./(sum(w).^2 - sum(w.^2)) * (w.*x)'*x; %weighted cov

