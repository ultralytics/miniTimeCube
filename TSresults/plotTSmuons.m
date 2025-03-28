% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = plotTSmuons()
close(findobj(0,'type','figure'))
if nargin==0
    [filename,pathname]=uigetfile([pwd '/TSresults/*.*'],'Select MC file to plot:'); if filename==0; fprintf('No file selected ... Done.\n'); return; end; addpath(pathname);
    fprintf('Loading MC file ''%s'' \nfrom folder ''%s''...',filename,pathname);  
    tic; load([pathname filename]); fprintf(' Done. (%.1fs)\n',toc)
end
MTCflag = false;


if ~isempty(strfind(MC.FileName,'Fiber'))
    i = ~all(MC.xhat==0,2) & ~any(isnan(MC.xhat),2);
    X = MC.xhat(i,:);  X = [fcnel(X(:,1:3))*r2d, fcnaz(X(:,1:3))*r2d, X(:,[4 5 1:3])]; %[el az de dx uvec];
    T = MC.xtrue(i,:); T = [fcnel(T(:,1:3))*r2d, fcnaz(T(:,1:3))*r2d, T(:,[4 5 1:3])]; %[el az de dx uvec];
    E = [X(:,1:4)-T(:,1:4) sum(X(:,5:7).*T(:,5:7),2)]; %[el az dx de ct]
    %E(:,2) = acosd(sum(fcnvec2uvec(MC.xtrue(i,1:2)).*fcnvec2uvec(MC.xhat(i,1:2)),2));
    
    ha=fig(2,4,'19x38cm'); nb=100; clear S;  s={'el (deg)','az (deg)','dE (MeV)','dx (mm)'};  ha(2).Title.String=str_(MC.FileName);
    j=X(:,3)>5 & X(:,4)>30 & X(:,4)<norm(input.cube.Lr*2*.92);  xa=X(j,1); std(E(j,:))
    mean(j)
    for i=1:4
        sca; hist2(T(j,i),X(j,i),nb);    xyzlabel(s{i},'fit')
        %sca; y=E(j,i);  hist2(xa,y,nb);   [~,~,S(i)]=movingMean(xa,y,nb,1);  xyzlabel(s{i},'fit')
    end
    for i=1:4
        sca; fcnhist(E(j,i),nb);    xyzlabel(s{i},'fit'); axis tight
    end
    fig; histogram(E(:,5),100);
        
    %fig(1,3); for i=1:3; sca(i); x=linspace(0,3,1000); plot(x,S(i).s(x),'.-','Display','CFD - 50%');  xyzlabel('E (MeV)',s{i}); end
    return
else
    MTCflag = true;
    dedx = 10*MC.xhat(:,3)./MC.xhat(:,4); %MeV/cm
    if MTCflag
        i = ~MC.failure(:,1) & MC.xhat(:,4)~=0 & dedx>0;
    else
        i = ~MC.failure(:,1) & MC.xtrue(:,3)~=0 & dedx>1.7 & dedx<3.5;
    end
    
    X = MC.xhat(i,:);
    T = MC.xtrue(i,:);
end


E = X-T;
dEdx = T(:,3)./T(:,4);
dEdxhat = X(:,3)./X(:,4);
if MTCflag
    E(:,5) = dEdxhat;
    dtheta = fcnangle(fcnSC2CC(1,X(:,1),X(:,2)),[0 0 -1]);
else
    E(:,5) = dEdxhat-dEdx;
    dtheta = fcnangle(fcnSC2CC(1,X(:,1),X(:,2)),fcnSC2CC(1,T(:,1),T(:,2)));
end
eff = sum(i)/numel(i);


n = 30;
if MTCflag
    fig(2,2,'20x20cm'); 
    sca;  fcnhistc(cos(dtheta).^2,n); xyzlabel('muon angle cos^2(\Theta)');
    sca;  a=E(:,3); fcnhistc(a,n); xyzlabel('muon dE (MeV)'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
    sca;  a=E(:,4); fcnhistc(a,n); xyzlabel('muon dx (cm)'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
    sca;  a=E(:,5); fcnhistc(a,n); xyzlabel('muon dE/dx (MeV/cm)'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
    fcnfontsize(16); fcntight('y')
else
    fig(2,4,1.1); 
    sca;  x=linspace(0,180,n); fcnhist(dtheta*r2d,x,'r');  %x=fcndsearch(dtheta*r2d,.683);
    %xyzlabel('angle (deg)','','',sprintf('''%s'' angles  (%.2f efficiency)',str_(input.MTC.runname),eff)); set(gca,'xlim',[-1 180]);
    xyzlabel('angle (deg)','','',sprintf('''%s'' angles  (%.2f efficiency)',str_(MC.FileName),eff)); set(gca,'xlim',[-1 180]);
    sca;  plot(T(:,3),X(:,3),'b.'); xyzlabel('true (MeV)','estimated (MeV)','','muon E');
    sca;  plot(T(:,4)*10,X(:,4)*10,'b.'); xyzlabel('true (cm)','estimated (cm)','','muon dx');
    sca;  plot(dEdx*10,dEdxhat*10,'b.'); xyzlabel('true (MeV/cm)','estimated (MeV/cm)','','muon dE/dx');
    sca;  fcnhistc(-cos(dtheta),n); xyzlabel('estimated - true cos(\Theta)','','','muon angle');
    sca;  a=E(:,3); fcnhistc(a,n); xyzlabel('estimated - true (MeV)','','','muon E'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
    sca;  a=E(:,4); fcnhistc(a,n); xyzlabel('estimated - true (cm)','','','muon dx'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
    sca;  a=E(:,5); fcnhistc(a,n); xyzlabel('estimated - true (MeV/cm)','','','muon dE/dx'); legend(sprintf('%.1f +/- %.1f',mean(a),std(a)))
end
