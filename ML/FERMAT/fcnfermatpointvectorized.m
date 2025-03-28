% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [fx, gx, F, sourceprobability, F1] = fcnfermatpointvectorized(k,X)
%X = [xyztw...xyztw]
gx      = []; %gradient
[nv,np] = size(X);  np=np/5; %nv = number of vectorized points
iv      = (0:np-1)'*5;
posi    = (1:3) + iv;

p = reshape(X(:,posi),[np*nv 3]);
[~, I, J] = fcnunique(p(:,1:3),'rows');

[Fs, r] = fcnsolidanglevectorized(k,p(I,:));
Fs = Fs(:,J,:).*fcnrow(X(:,iv+5)); %Fspatial

if k.cherenkov
    %r1 = r(:,:,1);
    [cr,cx] = rangec(X(:,posi(1,:)), X(:,posi(2,:)));  ux=cx./cr;
    ux = permute(repmat(ux,[1 1 np]),[1 3 2]);
    ux = reshape(ux,[np*nv 3]);

    %point travel uvectors
    [pr,~,px,py,pz] = bsxrangec(k.pxyz, p); %pixel-point vectors
    ct = (px.*ux(:,1)'+py.*ux(:,2)'+pz.*ux(:,3)')./pr;
    
    FCs = interp1c(k.smearExp.ax,k.smearExp.apdf,ct); %FCherenkovSpatial
end

k.nhpp(k.nhpp<.3)=0;
if k.timeflag
    nt = numel(k.ztpid);
    rt = k.zt(k.ztpid) - r(k.ztpid,:,:)*k.ci;
    dt = rt(:,J,:) - fcnrow(X(:,iv+4));
    Ft = interp1c(k.smearExp.x,k.smearExp.pdfsmeared,dt); %Ftemporal
    %F = sum(Ft.*Fs,3);
    
    if k.cherenkov
        Ft = interp1c(k.smearExp.cx,k.smearExp.cpdf,dt(:,:,1)); %FCherenkovTemporal
        Fs = Fs(:,:,1).*FCs;
        F = F*.96 + Ft.*Fs*.04;
    end
    %F1  = sum(reshape(F,[k.npixels nv np]),3); %sum sources

    Ft = sum(Ft,3); %sum reflections
    Fs = sum(Fs,3); %sum reflections
    
    Ft = reshape(Ft,[nt nv np]);
    Fs = reshape(Fs,[k.npixels nv np]);
    if np>1;  w=Fs(k.ztpid,:,:);  w=w./(sum(w,3)+1E-323);  Ft=Ft.*w;  end; %weigh Time likelihood by spatial ratio
    
    Ft  = sum(Ft,3); %sum sources
    lambda  = sum(Fs,3)*k.yield + 1E-3; %sum sources
    FP = exp( k.nhpp.*log(lambda) - (lambda+k.poiss.gammaln) ); %NLL from MATLAB poisspdf function. SLOWER, WORKS WITH NON INTEGER ZN's

    F=FP;  F(k.ztpid,:)=F(k.ztpid,:).*Ft;
    fx = -sum(log( max(F,1E-323) ));
else
    if k.cherenkov
        Fs = FCs; %Fs(:,:,1).*(.96 + FCs*.04);
    end
    F = sum(Fs,3); %sum reflections
    lambda  = sum(reshape(F*k.yield,[k.npixels nv np]),3) + 1E-3; %sum sources
    fx = -sum( k.nhpp.*log(lambda) - (lambda+k.poiss.gammaln) ); %NLL from MATLAB poisspdf function. SLOWER, WORKS WITH NON INTEGER ZN's

    %F1 = sum(reshape(F,[k.npixels nv np]), 3); %sum sources
    %fx = -sum(F1) %OLD WAY!!
end
%fx = -sum(F1);
%floor=1E-323*0;
%fx = -sum(log(F1+floor));
%fx = -k.nhpp'*log(F1+floor);    %fx  = -sumlogc(F2);

fx=fx(:);


if nargout>2
    sourceprobability = (F+1E-323)./sum(F+1E-323,2);
    if nargout==5;  F1=lambda; end
end


% 
% a = input.table.smearExp
% b = interp1c(a.cx,a.cpdf,a.x);
% fig; plot(a.x,a.pdfsmeared*.96+b*.04);