function neutronpdfplots(input,G1,handles,x,xcapture)
mnp = 3; %min number of photons
nb = 2; %size(x,3); %number of bounces
nt = size(x,1); %number of total samples
for b = 1:nb; %bounce
    i = x(:,1,b)==2212; %protons
    j = ~i & x(:,1,b)~=0; %C12, C13, etc
    inside = x(:,8,b)==1;
    dEhat = x(:,15,b)~=0; %dEhat>0
    observedPhotons = x(:,7,b)>mnp;
    %xp{b} = x(i & inside & dEhat & observedPhotons,:,b); %proton bounces
    %xc{b} = x(j & inside & dEhat & observedPhotons,:,b); %carbon bounces
    xp{b} = x(i,:,b); %proton bounces
    xc{b} = x(j,:,b); %carbon bounces
    
    xpm(b,:) = [b sum(i)/nt mean(x(i,:,b))]; %proton bounces
    xcm(b,:) = [b sum(j)/nt mean(x(j,:,b))]; %carbon bounces
end

%x = [1-pid, 2-neutron range (mm), 3-time (ns), 4-angle (rad), 5-dE (MeV), 6-emitted photons, 7-observed photons, 8-inside detector, 9-x, 10-y, 11-z, ...
%     12-xhat, 13-yhat, 14-zhat, 15-dEhat]
%xcapture = [neutroncapturexyz, 4-capturetime, 5-capturephotons, 6-allinsideflag, 7-E0, 8-E0enclosed, 9-E0hat]

clc
fprintf('Proton Recoils:\n')
fprintf('bounce    probability    GEANT pid   range (mm)    time (ns)  angle (deg)     dE (MeV)   emitted photons   collected photons  inside\n')
fprintf('%6.0f %14.2f %12.0f %12.1f %12.1f %12.1f %12.1f  %16.1f    %16.1f %7.2f\n',xpm(:,1:10)')
fprintf('\nCarbon Recoils:\n')
fprintf('bounce    probability    GEANT pid   range (mm)    time (ns)  angle (deg)     dE (MeV)   emitted photons   collected photons  inside\n')
fprintf('%6.0f %14.2f %12.0f %12.1f %12.1f %12.1f %12.1f  %16.1f    %16.1f %7.2f\n',xcm(:,1:10)')

fprintf('\n\nCandidate Criteria            Efficiency\n')
fprintf('%-30s%.2f\n','First Bounce on Proton',mean( x(:,1,1)==2212))
fprintf('%-30s%.2f\n','Second Bounce on Proton',mean( x(:,1,2)==2212 ))

v1 = x(:,1,1)==2212 & x(:,1,2)==2212;  fprintf('%-30s%.2f\n','P-P',mean(v1))
v2 = x(:,1,1)==2212 & x(:,1,2)~=2212; fprintf('%-30s%.2f\n','P-C',mean(v2))
v3 = x(:,1,1)~=2212 & x(:,1,2)==2212; fprintf('%-30s%.2f\n','C-P',mean(v3))
v4 = x(:,1,1)~=2212 & x(:,1,2)~=2212; fprintf('%-30s%.2f\n\n','C-C',mean(v4))

i = xcapture(:,6)==1;
fprintf('Candidate Criteria            P-P    P-C    C-P    C-C    All\n')
j=i;                        fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','Capture Inside Detector',mean(j(v1)),mean(j(v2)),mean(j(v3)),mean(j(v4)),mean(j))
j=i & xcapture(:,5)>2;      fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','Observe >2 Capture Photons',mean(j(v1)),mean(j(v2)),mean(j(v3)),mean(j(v4)),mean(j))
j=xcapture(:,4)<12000;      fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','Capture Time < 12us',mean(j(v1)),mean(j(v2)),mean(j(v3)),mean(j(v4)),mean(j))
j=x(:,5,1)>.100;             fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','dE01 > 100keV',mean(j(v1 & i)),mean(j(v2 & i)),mean(j(v3 & i)),mean(j(v4 & i)),mean(j(i)))
j=x(:,7,1)>10;              fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','bounce1 >10 photons',mean(j(v1 & i)),mean(j(v2 & i)),mean(j(v3 & i)),mean(j(v4 & i)),mean(j(i)))
j=x(:,7,2)>10;              fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','bounce2 >10 photons',mean(j(v1 & i)),mean(j(v2 & i)),mean(j(v3 & i)),mean(j(v4 & i)),mean(j(i)))
j=(x(:,5,1)./x(:,5,2))>.20;   fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','dE1/dE2 > 0.20',mean(j(v1 & i)),mean(j(v2 & i)),mean(j(v3 & i)),mean(j(v4 & i)),mean(j(i)))
j=xcapture(:,5)>2 & xcapture(:,6)==1 & xcapture(:,4)<12000 & x(:,5,1)>.100 & x(:,7,1)>10 & x(:,7,2)>10 & (x(:,5,1)./x(:,5,2))>.20;   
fprintf('%-30s%.2f%7.2f%7.2f%7.2f%7.2f\n','Combined',mean(j(v1)),mean(j(v2)),mean(j(v3)),mean(j(v4)),mean(j))


%PLOT STATS
close all
titles = {'Neutron Range (mm)','dt (ns)','d\Theta (deg)','dE (MeV)','Emitted Photons','Observed Photons (PE)'};
units = {'mm','ns','deg','MeV','photons','photons'};
bouncestr = {'First Bounce ','Second Bounce ','Third Bounce '};
for b = 1:nb; %bounce
    p = xp{b};
    c = xc{b};
    
    h=fig(4,2);
    [cstr, ccolor] = fcnpid2name(mode(c(:,1)));  ccolor = ccolor{1};
    if ~isempty(p)
        for i=1:6
            sca(h(i));
            [y,xv]=fcnhist(p(:,i+1),30); bar(xv,y,1,'edgecolor','b','facecolor',[.7 .7 1]);
            [y,xv]=fcnhist(c(:,i+1),30); bar(xv,y,1,'edgecolor',ccolor,'facecolor',ccolor + (1-ccolor)*.7);  
            axis tight; grid on; ylabel('Events'); xlabel(titles{i}); title([bouncestr{b} titles{i}]); alpha(0.8)
            
            legend(sprintf('Proton: %.1f+/-%.1f %s',mean(p(:,i+1)),std(p(:,i+1)),units{i}),sprintf('Carbon: %.1f+/-%.1f %s',mean(c(:,i+1)),std(c(:,i+1)),units{i}));
            if i==1; title(sprintf('%s\n%.1f%% Proton Recoils and %.1f%% Carbon Recoils',[bouncestr{b} titles{i}],xpm(b,2)*100,xcm(b,2)*100)); end
        end
    end
    
end


% %PLOT BOUNCE-ENERGY RELATIONSHIP
% newfitflag = 0;
% if ~newfitflag; load dEfit; end
% for b=1:nb
%     h=fig(4,2);  p = xp{b};  c = xc{b};
%     
%     sca(h(1))
%     plot(p(:,5),p(:,4),'b.'); hold on
%     plot(c(:,5),c(:,4),'.','color',ccolor); axis tight; grid on; ylabel('true d\Theta (deg)'); xlabel('true dE (MeV)'); title([bouncestr{b} 'true d\Theta vs true dE'])
%     legend('Proton','Carbon'); axis1 = axis;
%     
%     sca(h(3))
%     plot(p(:,6),p(:,4),'b.'); hold on
%     plot(c(:,6),c(:,4),'.','color',ccolor); axis tight; grid on; ylabel('true d\Theta (deg)'); xlabel('Emitted Photons'); title([bouncestr{b} 'true d\Theta vs Emitted Photons'])
%     
%     sca(h(5))
%     plot(p(:,7),p(:,4),'b.'); hold on
%     plot(c(:,7),c(:,4),'.','color',ccolor); axis tight; grid on; ylabel('true d\Theta (deg)'); xlabel('Observed Photons (PE)'); title([bouncestr{b} 'true d\Theta vs Observed Photons'])
%     
%     sca(h(2))
%     plot(p(:,15),p(:,5),'b.'); hold on
%     plot(c(:,15),c(:,5),'.','color',ccolor); ylabel('true dE (MeV)'); xlabel('estimated dE (MeV)'); title([bouncestr{b} 'true dE vs estimated dE'])
%     axis([0 max(p(:,15)) 0 max(p(:,5))]); grid on;
% 
%     if b==1 && newfitflag; dEfit = createFit(p(:,15), p(:,5)); end
%     x1=linspace(eps,max(p(:,15)),1000); plot(x1,dEfit(x1),'-','linewidth',2,'color',[.7 .7 .7])
%     e=[p(:,5),  (p(:,15)-p(:,5))./p(:,5),   (dEfit(p(:,15))-p(:,5))./p(:,5)];
% 
%     sca(h(6))
%     xv = linspace(-.8,.8,31);
%     bar(xv,hist(e(:,3),xv),1,'edgecolor',[.7 .7 1],'facecolor',[0 0 1]); hold on; ylabel('Events'); xlabel('error (fraction)'); title([bouncestr{b} 'energy estimate error (calibrated)'])
%     axis tight; grid on; legend(sprintf('%.2f 1\\sigma',fcnstdnonlin(e(:,3),[],median(e(:,3)))))
%     
%     sca(h(8))
%     e = p(:,12:14) - p(:,9:11);
%     xv = linspace(-25,25,31);
%     bar(xv,hist(e(:,1),xv),1,'edgecolor','r','facecolor',[1  .7 .7]); hold on;
%     bar(xv,hist(e(:,2),xv),1,'edgecolor','g','facecolor',[.7  1 .7]); hold on;
%     bar(xv,hist(e(:,3),xv),1,'edgecolor','b','facecolor',[.7 .7  1]); hold on;  ylabel('Events'); xlabel('error (mm)'); title([bouncestr{b} 'position estimate error'])
%     axis tight; grid on; legend(sprintf('x %.1fmm 1\\sigma',fcnstdnonlin(e(:,1))),sprintf('y %.1fmm 1\\sigma',fcnstdnonlin(e(:,2))),sprintf('z %.1fmm 1\\sigma',fcnstdnonlin(e(:,3))))
%     
%     sca(h(7))
%     plot(dEfit(p(:,15)),p(:,4),'.','color',[0 0 1]); hold on
%     plot(dEfit(c(:,15)),c(:,4),'.','color',ccolor); 
%     plot(p(:,15),p(:,4),'.','color',fcnlightcolor([0 0 1])); 
%     plot(c(:,15),c(:,4),'.','color',fcnlightcolor(ccolor));
%     axis tight; grid on; ylabel('true d\Theta (deg)'); xlabel('estimated dE (MeV)'); title([bouncestr{b} 'true d\Theta vs estimated dE'])
%     legend('Calibrated Proton','Calibrated Carbon','Proton','Carbon')
%     axis(axis1)
% end
% 
% 
% 
% %PLOT CAPTURE STUFF
% h = fig(2,2);
% 
% sca(h(1))
% a=xcapture(:,4); [y,xv]=fcnhist(a,50); bar(xv,y,1,'edgecolor',[.7 .7 1],'facecolor',[0 0 1]); axis tight; grid on; ylabel('Events'); xlabel('t (ns)'); title(sprintf('neutron capture time (%.3f < 12\\mus)',mean(a<12000))); hold on
% plot([1 1]*12E3,[0 max(y)],'b-','linewidth',2,'color',[.7 .7 .7])
% text(12E3,max(y)*.7,'12\mus','rotation',90,'fontsize',12,'verticalalignment','bottom')
% 
% sca(h(2))
% popoutsubplot(handles.GUI.axes1,h(2)); hold on
% set(findobj(gca,'Type','Line'),'LineWidth',2);
% a=xcapture(:,1:3); plot3(a(:,1),a(:,2),a(:,3),'b.','markersize',3); L=65; axis([-L L -L L -L L]*2.5)
% i= fcninsidedetector(a, input);
% title(sprintf('neutron capture pos (%.2f end inside, %.2f always inside)',mean(i),mean(i & xcapture(:,5)>2))); 
% set(gca,'XTickMode','auto','YTickMode','auto','ZTickMode','auto','XTickLabelMode','auto','YTickLabelMode','auto','ZTickLabelMode','auto')
% 
% sca(h(3)); cla
% i = xcapture(:,7)>.1 & xcapture(:,9)>.1 & x(:,1,1)==2212 & x(:,7,1)>10 & x(:,8,1) & x(:,1,2)==2212 & x(:,7,2)>10 & x(:,8,2);
% E0 = xcapture(i,7);
% E0e = xcapture(i,8); %true E0enclosed
% E0ehat = xcapture(i,9); %estimated E0enclosed
% if newfitflag; dEsumfit = createFit(E0ehat, E0e);  save('dEfit.mat','dEfit','dEsumfit'); end
% E0ehatc = dEsumfit(E0ehat); %calibrated
% plot(E0ehat,E0e,'.','color',fcnlightcolor([0 0 1])); hold on
% plot(E0ehatc,E0e,'b.'); axis tight
% plot(get(gca,'xlim'),[1 1]*E0(1),'-','linewidth',2,'color',[.7 .7 .7])
% axis tight; grid on; ylabel('E0e (MeV)'); xlabel('estimated E0e (MeV)'); title('E0 vs estimated E0');
% x1 = linspace(eps,max(E0ehat),1000);  plot(x1,dEsumfit(x1),'-','color','g','linewidth',1);
% legend('uncalibrated','calibrated','true E0','fit')
% 
% %save('x.mat','E0ehat','E0e')
% 
% sca(h(4)); cla
% a = (E0ehatc - E0)./E0;
% b = (E0ehatc - E0e)./E0;
% xv = linspace(-1,.8,31);
% bar(xv,hist(a,xv),1,'edgecolor',[1 .7 .7],'facecolor',[1 0 0]); hold on; ylabel('Events'); xlabel('error (fraction)'); title('E0 = sum(dE) estimate error')
% bar(xv,hist(b,xv),1,'edgecolor',[.7 .7 1],'facecolor',[0 0 1]);
% axis tight; grid on; legend(sprintf('E0 error, %.2f 1\\sigma',fcnstdnonlin(a,[],0)),sprintf('E0 enclosed error, %.2f 1\\sigma',fcnstdnonlin(b,[],0)))


% %CORRELATIONS
% %[1-pid, 2-neutron range (mm), 3-time (ns), 4-angle (rad), 5-dE (MeV), 6-emitted photons, 7-observed photons, 8-inside detector, 9-x, 10-y, 11-z, ...
% %     12-xhat, 13-yhat, 14-zhat, 15-dEhat, 16-dealyed bookend observed photons]
% 
% %1-range, 2-time, 3-angle, 4-dE, 5-emitted photons, 6-observed photons, 7-estimated dE
% mnp = 10;
% i = x(:,1,1)==2212 & x(:,1,2)==2212; %protons
% inside = x(:,8,1)==1 & x(:,8,2)==1;
% converged = x(:,15,1)~=0 & x(:,7,1)>mnp & x(:,15,2)~=0 & x(:,7,2)>mnp;
% a = x(i & inside & converged,:,1);  a = a(:,[2:7 15]);  %a(:,3) = sin(a(:,3)/57.3).^2;
% b = x(i & inside & converged,:,2);  b = b(:,[2:7 15]);  %b(:,3) = sin(b(:,3)/57.3).^2;
% titles{6} = sprintf('Observed Photons');
% titles{7} = sprintf('Estimated\n\dE (MeV)');
% 
% fig;  [~,h,~,~,p] = plotmatrix(a); title('First Proton Bounce Correlations')
% for i=1:size(a,2)
%    set(get(h(end,i),'XLabel'),'String',titles{i})
%    set(get(h(i,1),'YLabel'),'String',titles{i})
% end; axis(h(:),'tight'); axis(p,'tight'); fcnmarkersize(3)
% 
% 
% fig;  [~,h,~,~,p] = plotmatrix(b); title('Second Proton Bounce Correlations')
% for i=1:size(a,2)
%    set(get(h(end,i),'XLabel'),'String',titles{i})
%    set(get(h(i,1),'YLabel'),'String',titles{i})
% end; axis(h(:),'tight'); axis(p,'tight'); fcnmarkersize(3)
% 
% fig;  [~,h,~,~,p] = plotmatrix(a,b); title('First-Second Proton Bounce Correlations')
% for i=1:size(a,2)
%    set(get(h(end,i),'XLabel'),'String',titles{i})
%    set(get(h(i,1),'YLabel'),'String',titles{i})
% end; axis(h(:),'tight'); axis(p,'tight'); fcnmarkersize(3)
% 
% fig; [~,h,~,~,p] = plotmatrix([a b]);
% axis(h(:),'tight'); axis(p,'tight')
% h1 = h(1:size(a,2),size(a,2)+1:size(a,2)*2); h1=cell2mat(get(h1(:),'children')); set(h1,'color',[.7 .7 .7])
% h1 = h(size(a,2)+1:size(a,2)*2,1:size(a,2)); h1=cell2mat(get(h1(:),'children')); set(h1,'color',[.7 .7 .7])
% title('First-Second Proton Bounces, Full Cross Correlation Matrix'); fcnmarkersize(2)
end

function c = fcnlightcolor(c)
c = (1 - c)*.7 + c;
end


function [fitresult, gof] = createFit(x1, y1)
[xData, yData] = prepareCurveData( x1, y1 );

% % Set up fittype and options.
% ft = fittype( 'power1' );
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.Lower = [-Inf -Inf];
% opts.Robust = 'LAR';
% opts.StartPoint = [223.524193560009 1.66725570194361];
% opts.Upper = [Inf Inf];

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf];
opts.Robust = 'LAR';
opts.StartPoint = [223.524193560009 1.66725570194361 -0.684419446848861];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y1 vs. x1', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'x1' );
% ylabel( 'y1' );
% grid on
end

