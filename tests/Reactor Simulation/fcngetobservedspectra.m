function s = fcngetobservedspectra(shield,plotflag)
if nargin==1; plotflag=1; end
close all
colors = {[1 0 0],[0 1 0],[0 0 1],[1 .7 0],[1 0 1],[.8 .8 .8],[0 1 1],[.5 0 1]};
if nargin == 0
    shield.lead = 0;
    shield.polyethylene = 0;
end

s = fcnenergydistributions(); %spectra
%c = getcandidatecuts(); %candidates

fnames = {'MC-2L PLANACON XP85012-A1 80pc 1536pix 29qe 13re 2cp 5%EJ-254 10k Antineutrino1000',....
    'MC-2L PLANACON XP85012-A1 80pc 1536pix 29qe 13re 2cp 5%EJ-254 10k Neutron',...
    'MC-2L PLANACON XP85012-A1 80pc 1536pix 29qe 13re 2cp 5%EJ-254 10k Gamma',...
    'MC-2L PLANACON XP85012-A1 80pc 1536pix 29qe 13re 2cp 5%EJ-254 10k Muon'};
fi=0;
for i=1:6
    if any(i==[1 2 4 6])
        fi=fi+1;
        MC = getcandidatefractions(fnames{fi});
        xstr = sprintf('%s Energy (MeV)',MC.particleName);
    end
    legendstr = {[MC.particleName ' Spectrum'],'Full Event','Prompt -> Prompt','Prompt -> Delayed','Delayed -> Prompt','Delayed -> Delayed','Unobserved','Rejected','Shielded'};
    
    
    switch i
        case 1 %antineutrinos
            MC.fraction(:,end+1) = 0;
        case 2 %reactor neutrons
            MC.fraction(:,end+1) = 0;
        case 3 %atmospheric neutrons
            MC.fraction(:,end+1) = 0;
        case {4,5} %reactor and atmospheric gammas
            MC = getcandidatefractions(fnames{fi});
            
            pf = getGammaTransInLead(MC.MeVVector, shield.lead)'; sf = 1-pf;  %pass and shielded fractions (5cm lead)
            MC.fraction = MC.fraction.*pf;
            MC.fraction(:,end+1) = sf;
            legendstr{9} = sprintf('Shielded %.0fcm lead',shield.lead);
        case 6 %muons
            MC.fraction(:,end+1) = 0;
    end
    
    sx=s(i).pdfx;  sy=s(i).pdfy;  sy=sy/sum(sy)*100;
    
    s(i).pdfyc = max(interp1(MC.MeVVector',MC.fraction,sx,'linear',0), 0); %candidates
    s(i).sf = sy*s(i).pdfyc; %scalar fraction
    s(i).simplename = s(i).name(1:find(s(i).name==':')-1);
    
    if plotflag
        h=fig(2,2,'28.5cm');
        sca(h(1));
        loglog(sx,sy,'k','linewidth',2); hold on;
        for j=1:8
            legendstr{j+1} = sprintf('%s %.2g%%',legendstr{j+1},s(i).sf(j));
            plot(sx,s(i).pdfyc(:,j)','color',colors{j},'linewidth',1.5);
        end
        xlim = sx([find(sy>0,1,'first') find(sy>0,1,'last')]);
        grid on; axis tight; set(gca,'xlim',xlim); xlabel(xstr); ylabel('Candidate Fraction')
        %title(sprintf('%s Candidate Fractions',MC.particleName))
        title(sprintf('%s',s(i).name))
        legend(legendstr,'Location','SouthWest')
        
        
        sca(h(3))
        loglog(sx,sy,'k','linewidth',2); hold on;
        for j=1:8
            plot(sx,s(i).pdfyc(:,j)'.*sy,'color',colors{j},'linewidth',1.5);
        end
        grid on; axis tight; set(gca,'xlim',xlim); xlabel(xstr); ylabel('pdf')
        title(sprintf('%s Candidate Spectra',MC.particleName))
        legend(legendstr,'Location','SouthWest')
        
        
        sca(h(2))
        area(sx',s(i).pdfyc(:,1:8).*(sy'*ones(1,8)))
        
        ch = get(gca,'Children');
        for j=1:8
            set(ch(9-j),'FaceColor',colors{j},'edgecolor','none')
        end
        grid on; axis tight; set(gca,'xlim',xlim); xlabel(xstr); ylabel('pdf')
        set(gca,'xscale','log')
        title(sprintf('%s Complete Spectrum',MC.particleName))
        legend(legendstr(2:9),'Location','NorthEast')
        
        
        sca(h(4)); axis off; daspect([1 1 1])
        y=s(i).sf(1:8)/100; j=find(y>.001);  nj=numel(j);  vls=2:9;
        if nj>0
            pie(y(j),ones(1,nj),legendstr(vls(j)))
            ch = get(gca,'Children');
            l=nj*2+2;
            for k=j
                l=l-2;
                set(ch(l),'FaceColor',colors{k},'edgecolor','none')
                set(ch(l-1),'color',[.5 .5 .5],'Margin',2)
            end
        end
        set(h(1:3),'yscale','log','xscale','log')
        fcnfontsize(12)
    end
end
end

function MC = getcandidatefractions(fname)
input=[];
load(fname)
ni = length(MeVVector);  MC.MeVVector=MeVVector;
nj = numel(MC.xhat(:,1,1)); %#ok<*NODEF>

zv  = zeros(nj,ni);
combinedPrompt = zv;
combinedDelayed = zv;
combinedEvent = zv;
for pdi=1:2; %prompt-delayed index. 1=prompt, 2=delayed
    cp = false(nj,ni);
    cd = false(nj,ni);
    for i=1:ni
        pwc = 5; %(mm) positron wall distance [min]
        nwc = 5; %(mm) neutron wall distance [min]
        pndc = 250; %(mm) positron-neutron distance [max]
        pec = [1.000 8.000]; %(MeV) prompt energy [min max]
        dec = [0.050 0.100]; %(MeV) delayed energy [min max]
        dtc = [0.050 12.000]; %(us) delay time [min max]
        Ep = MC.xhat(:,4,i);
        En = MC.xhat(:,8,i);
        xn = MC.xhat(:,5:7,i);
        xp = MC.xhat(:,1:3,i);
        dt = (MC.firstHit(:,2,i)-MC.firstHit(:,1,i))/1000;
        
        if pdi==1; %assume prompt signal is faking
            En=Ep;
            xn=xp;
        elseif pdi==3 %assume delayed signal is faking
            Ep = En;
            xp = xn;
        end
        
        positron_energy1 = Ep>pec(1) & Ep<pec(2);
        neutron_energy1 =  En>dec(1) & En<dec(2);
        positron_geode_distance1 = fcninsidedetector(xp, input, pwc);
        neutron_geode_distance1 = fcninsidedetector(xn, input, nwc);
        neutron_delay1 = dt>dtc(1) & dt<dtc(2);
        positron_neutron_distance1 = rangec(xn, xp) < pndc;
        if strcmp(input.cube.shape,'cylinder')
            r1 = fcnrange(xp(:,1:2));
            r2 = fcnrange(xn(:,1:2));
            positron_geode_distance1 = abs(xp(:,3)) < (input.cube.height/2-pwc) & r1 < (input.cube.radius-pwc);
            neutron_geode_distance1 = abs(xn(:,3)) < (input.cube.height/2-nwc) & r2 < (input.cube.radius-nwc);
        end
        
        pc1 = ~MC.failedConvergenceFlag(:,pdi,i); %prompt convergence

        cp1 = (pc1.*positron_energy1.*positron_geode_distance1) ==1; %converged prompt
        cd1 = (pc1.*neutron_energy1.*neutron_geode_distance1) ==1; %converged delayed
        
        combinedPrompt(:,i,pdi) = cp1;
        combinedDelayed(:,i,pdi) = cd1;
        if pdi==2;  %full event
            combinedEvent(:,i) = neutron_delay1 & positron_neutron_distance1 & ~MC.failedConvergenceFlag(:,1,i) & pc1;
        end
        
        y1=double(MC.collectedPhotonCount(cp1,pdi,i));  plots.photons.mu(i,1) = mean(y1, 1);  plots.photons.s(i,1) = std(y1, 1);
        y2=double(MC.collectedPhotonCount(cd1,pdi,i));  plots.photons.mu(i,2) = mean(y2, 1);  plots.photons.s(i,2) = std(y2, 1);
        
        y1=Ep(cp1);  plots.Ehat.mu(i,1) = mean(y1, 1);  plots.Ehat.s(i,1) = std(y1, 1);
        y2=En(cd1);  plots.Ehat.mu(i,2) = mean(y2, 1);  plots.Ehat.s(i,2) = std(y2, 1);
        
        y1=MC.lastHit(cp1,pdi,i)-MC.firstHit(cp1,pdi,i);  plots.duration.mu(i,1) = mean(y1);  plots.duration.s(i,1) = std(y1);
        y2=MC.lastHit(cd1,pdi,i)-MC.firstHit(cd1,pdi,i);  plots.duration.mu(i,2) = mean(y2);  plots.duration.s(i,2) = std(y2);
        
        cp(:,i)=cp1;
        cd(:,i)=cd1;
    end
    MC.combinedPrompt(:,pdi) = sum(combinedPrompt(:,:,pdi))/nj;
    MC.combinedDelayed(:,pdi) = sum(combinedDelayed(:,:,pdi))/nj;
end
MC.combinedEvent = combinedPrompt(:,:,1) & combinedDelayed(:,:,2) & combinedEvent;

fraction = zeros(ni,7);
fraction(:,1) = sum(MC.combinedEvent)/nj; %full event
fraction(:,2) = sum(~MC.combinedEvent & combinedPrompt(:,:,1))/nj; %prompt faking prompt
fraction(:,3) = sum(~MC.combinedEvent & combinedDelayed(:,:,1))/nj; %prompt faking delayed
fraction(:,4) = sum(~MC.combinedEvent & combinedPrompt(:,:,2))/nj; %delayed faking prompt
fraction(:,5) = sum(~MC.combinedEvent & combinedDelayed(:,:,2))/nj; %delayed faking delayed
fraction(:,6) = reshape(sum(MC.collectedPhotonCount(:,1,:)==0 & MC.collectedPhotonCount(:,2,:)==0, 1)/nj,[ni 1]); %missed entire event
fraction(:,7) = ones(ni,1) - sum(fraction(:,1:6),2); %rejected entire event
%fraction(:,8) = sum(fraction(:,1:5),2); %some kind of candidate
MC.fraction=fraction;

%save([fcnfile2folder(MC.FileName) MC.FileName],'MC','MeVVector','nj','ni','input','vi','MCtable');
end