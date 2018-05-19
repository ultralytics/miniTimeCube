function [] = fcnPlotDendrogram(input,output,photons,GEANT,G1)
X=G1.dendrogram;
tid = G1.tid;
pid = G1.pid;
ptid = G1.ptid;
maxgen = size(X,2)-1;
scale = 1*max(sqrt(maxgen),2);
h = fig(1,1,scale,scale,[2, 2, 0, 0, 0.1, 1.5]);  fcn3label;  set(h,'ydir','reverse','zdir','reverse');  fcnview('top'); axis off

P = zeros(G1.nutid,3); %position of each point
I = zeros(G1.nutid,1); %index of parent
A = ones(G1.nutid,1)*2*pi; %arc angle available to particle
AZ = zeros(G1.nutid,1); %azimuth angle of particle
PID = fcntid2pid(X(:,1),tid,pid);
[name, color, mass, charge] = fcnpid2name(PID);  name=str_(name);

pnf = 1;  % process name flag
if isfield(G1,'processname')
    processname = fcntid2processname(X(:,1),G1);
else
    pnf = 0;
end


for i = 1:G1.nutid
    try %#ok<*TRYNC>
        I(i) = find(G1.utid==fcntid2ptid(G1.utid(i),tid,ptid)); %indices to parent row
    end
end


zv = zeros(G1.nutid,1);
htn=text(zv,zv,zv,name,'VerticalAlignment','Bottom','HorizontalAlignment','Right','fontsize',12); %names
if pnf; htp=text(zv,zv,zv,processname,'VerticalAlignment','Top','HorizontalAlignment','Right','fontsize',10); end %processes

r=0;
for i = 1:maxgen
    r = r+1; %radius
    v1 = find(X(:,i+1)==0);  n1=numel(v1);
    
    if i==1 %get smaller angles
        A(v1) = A(v1)/n1;  if n1==1; A(v1)=2*pi*.66; end
        az = linspace(0,2*pi,n1+1);  AZ(v1) = az(1:n1)+pi/2;
        P(v1,:) = fcnSC2CC(r,0,AZ(v1));
        
        for k=1:n1
            l = v1(k);
            x = P(l,:);
            plot3([x(1) 0],[x(2) 0],[x(3) 0],'-','color',color{v1(k)},'linewidth',3);
            
            ke1 = max(G1.ke1(tid==X(l,1)));
            str = sprintf('%s %s',dE2strdE(ke1),name{l});
            set(htn(l),'string',str,'position',x*.9,'color',color{l}*.8,'Rotation',90-AZ(l)*r2d)
            if pnf; set(htp(l),'position',x*.9,'color',color{l}*.8,'Rotation',90-AZ(l)*r2d); end
        end
    else
        upi = fcnunique(I(v1));  nupi=numel(upi); %unique parent indices
        
        for j=upi'
            v2 = find(I(v1)==j);  n2=numel(v2);  v1j = v1(v2);  sc=zeros(n2,3);  sc(:,1)=1;
            A(v1j) = A(j)/n2;
            
            
            da = A(j);
            az = (  (-da/2+da/(n2*2)) : da/n2 : (da/2-da/(n2*2))  )*1.3;
            %az = linspace(-A(j)/2,A(j)/2,n2);
            sc(:,3)=az;  cc=fcnSC2CC(sc);
            
            
            C = fcnRPY2DCM_B2W([0 0 AZ(j)]);
            cc_r = cc*C';  sc_r = fcnCC2SCr(cc_r); AZ(v1j)=sc_r(:,3);
            cc_w = cc_r + P(j,:);
            P(v1j,:) = cc_w;
            
            
            %pid(v1j)
            
            for k=1:n2
                l = v1j(k);
                x = [cc_w(k,:)' P(j,:)'];
                                
                
                if i<4 && n2<15
                    ke1 = max(G1.ke1(tid==X(l,1)));
                    str = sprintf('%s %s',dE2strdE(ke1),name{l});
                    set(htn(l),'string',str,'position',x(:,1)*.9+x(:,2)*.1,'color',color{l}*.8,'Rotation',90-AZ(l)*r2d)
                    if pnf; set(htp(l),'position',x(:,1)*.9+x(:,2)*.1,'color',color{l}*.8,'Rotation',90-AZ(l)*r2d); end
                else
                    set(htn(l),'position',x(:,1)*.9+x(:,2)*.1,'color',color{l}*.8,'Rotation',90-AZ(l)*r2d)
                    if pnf; deleteh(htp(l)); end
                end
            end
        end
        
    end
    
    if i<3
        hc = fcncircle(r,[0 0 0]); set(hc,'Color',[1 1 1]*.9,'marker','none')
    end
end


for i = 1:G1.nupid
    j = find(PID == G1.upid(i));  J=I(j);
    plot3(P(j,1),P(j,2),P(j,3),'.','color',G1.upidcolor{i},'markersize',25);
    
    k = J>0;
    if ~all(k); J=J(k); j=j(k); end
    plot3([P(J,1), P(j,1)]',[P(J,2), P(j,2)]',[P(J,3), P(j,3)]','-','color',G1.upidcolor{i})
end

fname = str_(GEANT.filename);
title(sprintf('GEANT Tree - event %g\n%s',input.eventNumber,fname),'fontsize',13)
fcntight 
set(gca,'DataAspectRatio',[1 1 1]); %axis equal
%fcnfontsize(10)

end


function processname = fcntid2processname(tid1,G1)
n = numel(tid1);
processname = cell(n,1); 

for i=1:n
    pn1 = G1.processname(G1.tid==tid1(i));
    [upn1, ~, B]= unique(pn1,'stable');
    nunique = accumarray(B,1,size(upn1))';
    for j = find(nunique>1)
        upn1{j} = sprintf('%s %gX',upn1{j},nunique(j));
    end

    processname{i} = sprintf('%s\n',upn1{:});
end

end


function str = dE2strdE(E)
%E in MeV in
if E>=1E6
    units = 'TeV';  E=E/1E6;    str = sprintf('%.3g %s',E,units);
elseif E>=1E3
    units = 'GeV';  E=E/1E3;    str = sprintf('%.3g %s',E,units);
elseif E>=1
    units = 'MeV';              str = sprintf('%.3g %s',E,units);
elseif E>=1E-3
    units = 'keV';  E=E*1E3;    str = sprintf('%3.0f %s',E,units);
elseif E>=1E-6
    units = 'eV';  E=E*1E6;    str = sprintf('%3.0f %s',E,units);
else
    str = sprintf(' <1 eV'); 
end
end
