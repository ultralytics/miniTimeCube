function h = fcnPlotPhotonEmissionTimes(input,output,photons,PE,handles,G1)
MTCflag = ischecked(handles.GUI.realdataflag);
if MTCflag
    h=fig(1,1,1.5);
    fcnhist(output(1).t,50,'b','bar')
    xyzlabel('t (ns)','','',str_(input.MTC.A.filename))
    return
end



h = fig(2,2,1.5);
[names, c] = fcnpid2name(G1.upid);
nb = 100; %number bins
promptPhotons = photons.sourceTime<25; %<25ns
promtpPE = PE.t<25;

v1=promptPhotons;
if sum(v1)>0
    t1=0;  
    t2=max3([photons.endTime(v1) photons.endTime(v1)])+.1;  
    x = linspace(t1,t2,nb);
    
    %PROMPT PHOTONS
    sca(h(1))
    j = zeros(1,G1.nupid);
    for i=1:G1.nupid
        t = photons.sourceTime(v1 & photons.Gppid==G1.upid(i));
        if ~isempty(t)
            j(i) = numel(t);
            n=hist(t,x);  %v1=n~=0;  n(v1)=log(n(v1));
            hb = bar(x,n,1,'FaceColor',c{i},'EdgeColor',c{i});  hold on;  set(get(hb,'Children'),'FaceAlpha',.5)
        end
    end
    nj = sum(j>0);  upid = G1.upid(j>0);  j=j(j>0);  lstr = cell(nj,1);  for i=1:nj; lstr{i} = sprintf('%.0f %s',j(i),names{i}); end
    axis tight; set(gca,'XLim',[t1 t2]); title(sprintf('%.0f Emitted Photons (prompt)',sum(j))); ylabel('photons'); xlabel('t (ns)'); legend(lstr)
    
    %PROMPT PE
    sca(h(3))
    j = zeros(1,G1.nupid);
    for i=find(G1.upidflag)
        t = PE.t(PE.Gpid==G1.upid(i) & promtpPE);
        if ~isempty(t)
            j(i) = numel(t);
            n=hist(t,x);  %v1=n~=0;  n(v1)=log(n(v1));
            hb = bar(x,n,1,'FaceColor',c{i},'EdgeColor',c{i});  hold on;  set(get(hb,'Children'),'FaceAlpha',.5)
        end
    end
    nj = sum(j>0);  upid = G1.upid(j>0);  j=j(j>0);  lstr = cell(nj,1);  for i=1:nj; lstr{i} = sprintf('%.0f %s',j(i),names{i}); end
    axis tight; set(gca,'XLim',[t1 t2]); title(sprintf('%.0f Observed Photons (prompt)',sum(j))); ylabel('photons'); xlabel('t (ns)');  legend(lstr)
else
    title(h(2),sprintf('%.0f Emitted Photons (delayed)',0)); ylabel('photons'); xlabel('t (ns)');
    title(h(4),sprintf('%.0f Emitted Photons (delayed)',0)); ylabel('photons'); xlabel('t (ns)');
end


v1=~promptPhotons;
if sum(v1)>0
    t1=min(photons.sourceTime(v1))-.1;  t2=max(photons.endTime(v1))+.1;  x = linspace(t1,t2,nb);
    
    %DELAYED PHOTONS
    sca(h(2))
    j = zeros(1,G1.nupid);
    for i=find(G1.upidflag)
        t = photons.sourceTime(v1 & photons.Gppid==G1.upid(i));
        if ~isempty(t)
            j(i) = numel(t);
            n=hist(t,x);  %v1=n~=0;  n(v1)=log(n(v1));
            hb = bar(x,n,1,'FaceColor',c{i},'EdgeColor',c{i});  hold on;  set(get(hb,'Children'),'FaceAlpha',.5)
        end
    end
    lstr = fcnpid2name(G1.upid(j>0));  j=j(j>0);  for i=1:numel(lstr); lstr{i} = sprintf('%.0f %s',j(i),lstr{i}); end
    axis tight; set(gca,'XLim',[t1 t2]); title(sprintf('%.0f Emitted Photons (delayed)',sum(j))); ylabel('photons'); xlabel('t (ns)'); legend(lstr)
    
    %DELAYED PE
    sca(h(4))
    j = zeros(1,G1.nupid);
    for i=find(G1.upidflag)
        t = PE.t(PE.Gpid==G1.upid(i) & ~promtpPE);
        if ~isempty(t)
            j(i) = numel(t);
            n=hist(t,x);  %v1=n~=0;  n(v1)=log(n(v1));
            hb = bar(x,n,1,'FaceColor',c{i},'EdgeColor',c{i});  hold on;  set(get(hb,'Children'),'FaceAlpha',.5)
        end
    end
    lstr = fcnpid2name(G1.upid(j>0));  j=j(j>0);  for i=1:numel(lstr); lstr{i} = sprintf('%.0f %s',j(i),lstr{i}); end
    axis tight; set(gca,'XLim',[t1 t2]); title(sprintf('%.0f Observed Photons (delayed)',sum(j))); ylabel('photons'); xlabel('t (ns)'); legend(lstr)
else
    title(h(2),sprintf('%.0f Emitted Photons (delayed)',0)); ylabel('photons'); xlabel('t (ns)');
    title(h(4),sprintf('%.0f Emitted Photons (delayed)',0)); ylabel('photons'); xlabel('t (ns)');
end
