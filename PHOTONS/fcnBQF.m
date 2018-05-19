function devis = fcnBQF(de,dx,BQF)
%converts energy to visible energy using BQF 
%dx = distance travelled in mm
%BQF = 0.126 mm/MeV
devis = de./(1+BQF*de./dx);



% %perform fine BQF integration
% devis1 = fcnBQF(de,r,BQF);
% n=9000;
% dt = (t2-t1)/(n-1);
% de2 = (ke1-ke2)/(n-1);
% ke = linspace(ke1,ke2,n);
% s = fcnke2c(ke(1:end-1),mass);
% r2 = s.*dt*299.792458; %c=300mm/ns
% devis2 = sum( fcnBQF(de2,r2*(r/sum(r2)),BQF)  );
% [devis1/devis2]

end

function plotGEANTdata()

%PLOT GEANT QUENCHINGS
G=GEANT;  n=size(G.x,1); i=1:n-1; j=2:n;  x=double(G.x);
badIndices=G.ei(1:end-1,2); 
j(badIndices)=[];
i(badIndices)=[];

event=x(j,1);
pos=x(:,5:7);
p1  = pos(i,:);
p2  = pos(j,:);
t1  = x(i,8);
t2  = x(j,8);
ke1 = x(i,9);
ke2 = x(j,9);
de  = x(j,10);
pid = x(j,2);
r=rangec(p1,p2);
[upid, ~, ~, nupid] = fcnunique(pid);

k=de>0 & r>0; ha=fig;
for ii=1:numel(upid)
    l=k & pid==upid(ii) & event==4;
    [name, color, mass, charge] = fcnpid2name(upid(ii));
    plot(de(l),r(l),'.','markersize',1,'Color',color{1},'DisplayName',name{1})
end
ha.XLim=[.01 3]; ha.YLim=[1 100];  xyzlabel('dE (MeV)','dx (mm)','',str_(G.filename)); fcnfontsize(14)
ha.YScale='log'; ha.XScale='log'; fcntight('xy')
end