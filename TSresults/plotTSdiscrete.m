function [] = plotTS(input,MC,tsv)
close(findobj(0,'type','figure'))
if nargin==0
    [filename,pathname]=uigetfile([pwd '/*.*'],'Select MC file to plot:'); if filename==0; fprintf('No file selected ... Done.\n'); return; end; addpath(pathname);
    fprintf('Loading MC file ''%s'' \nfrom folder ''%s''...',filename,pathname);  
    tic; load([pathname filename]); fprintf(' Done. (%.1fs)\n',toc)
end

X=MC.xhat;
nmc = size(X,1); %number of mc runs
nxp = size(X,2); %number of xhat parameters
nts = size(X,3); %number of ts runs

%xlstr = 'E_\nu (MeV)'; %x label string
%xlstr = 'MCP Count';
xlstr = 'Detector Side Length (m)';
%xlstr = 'Detector Pixels per Row';  tsv = tsv.^2*6;
%xhat = [x y z e t]

[C, XC, TC, EC, EFC, Cstr] = fcngetcandidates(input,MC.collectedPhotonCount,X,MC.xtrue,MC.failure,'antineutrino');


%PLOTS --------------------------------------------------------------------
h = fig(3,4,1,.75);
[~,colors] = fcndefaultcolors(1:3,3);

sca(h(1));  delete(h(5));  cp = get(gca,'Position'); set(gca,'Position',cp.*[1 .53 1 2.5])
for i=1:11
    plot(tsv,squeeze(mean(C(:,i,:),1)),'.-','Color',fcndefaultcolors(i,12));
end
ylabel('fraction'); title('Candidates')
hl=legend(Cstr,'Location','Best'); legend boxoff; set(hl,'fontsize',7)
axis([min(tsv) max(tsv)+.001 0 1])

sca(h(2))
for i = 1:3
    errorbar(tsv,EC.mu(:,i),EC.sigma(:,i),'.-','Color',colors{i});
end
ylabel('error (mm)');  title('Positron Start');  legend('X','Y','Z')

sca(h(3))
for i = 1:3
    errorbar(tsv,EC.mu(:,i+5),EC.sigma(:,i+5),'.-','Color',colors{i});
end
ylabel('error (mm)');  title('Neutron Capture');  legend('X','Y','Z')

sca(h(4))
volume = 1;%tsv'.^3;
efficiency = squeeze(mean(C(:,11,:),1));
plot(tsv,mean(EC.sigma(:,(1:3)+0),2).^-2.*volume.*efficiency,'.-','Color',colors{1});
plot(tsv,mean(EC.sigma(:,(1:3)+5),2).^-2.*volume.*efficiency,'.-','Color',colors{2});
set(gca,'yscale','log');
ylabel('mm^{-2}');  title('vertex information');  legend('prompt','delayed','Location','Best')

sca(h(6))
for i=1:nts
   j = C(:,11,i);
   vechat = X(j,(1:3)+5,i) - X(j,(1:3),i);
   vh.mu(i,1:3) = mean(vechat,1);
   vh.sigma(i,1:3) = std(vechat,[],1);
end
for i = 1:3
    errorbar(tsv,vh.mu(:,i),vh.sigma(:,i),'.-','Color',colors{i});
end
ylabel('error (mm)');  title('Direction Vector');  legend('X','Y','Z')

sca(h(7))
magnitude = fcnrange(vh.mu); %magnitude
sigma = mean(vh.sigma,2);
SNR = magnitude./sigma;
h1=bar(tsv,SNR,1,'b'); set(h1,'edgecolor',[.7 .7 .7]); 
ylabel('SNR'); title(sprintf('Direction Vector SNR (\\mu/\\sigma), %.4g mean',mean(SNR)))

sca(h(8))
plot(tsv,SNR.^2.*volume.*efficiency,'.-','Color',colors{3});
set(gca,'yscale','log');
ylabel('SNR^{-2}');  title('angle information');

sca(h(9))
z = double(MC.collectedPhotonCount);
mu = squeeze(mean(z,1))';  %sigma = squeeze(std(z,[],1))';
plot(tsv,mu(:,1),'.-','Color',colors{1});
plot(tsv,mu(:,2),'.-','Color',colors{2});  set(gca,'yscale','log');
%errorbar(tsv,mu(:,1),sigma(:,1),'.-','Color',colors{1},'MarkerSize',10);
%errorbar(tsv,mu(:,2),sigma(:,2),'.-','Color',colors{2},'MarkerSize',10);
ylabel('#');  title('observed photons');  legend('prompt','delayed','Location','Best')

sca(h(10))
mu = EC.mu(:,4);
sigma = EC.sigma(:,4);
errorbar(tsv,mu,sigma,'.-','Color',colors{1});
ylabel('error (MeV)'); title('prompt visible energy')
legend('MeV','Location','Best'); legend boxoff

sca(h(11))
i=isfinite(EC.mu(:,4));  tsvi=tsv(i);
%mu = EC.sigma(i,4);
%mu = mu'./mean(TC.mu(i,4));
mu = EFC.mu(i,4);
plot(tsvi,mu,'r.','MarkerSize',12,'LineWidth',1);
axis([min(tsvi) max(tsvi) 0 max(mu)*1.5]);
ylabel('error 1\sigma (fraction)');  title('prompt visible energy')

sca(h(12))
plot(tsv,EC.sigma(:,4+0).^-2.*volume.*efficiency,'.-','Color',colors{1});
ylabel('MeV^{-2}');  title('energy information');  legend('prompt','Location','Best')

%FIGURE-WIDE OPERATIONS
for i=1:numel(h)
   if isgraphics(h(i)); sca(h(i)); axis tight; end
end
%fcntight(h)
for i=[1 8 11 12]
   sca(h(i)); set(gca,'ylim',get(gca,'ylim').*[0 1]);
end






end

