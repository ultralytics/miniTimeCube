function [photons, PE]= fcnphotontransport(input, P, G1, flags)
%fprintf('Photons... '); startclock=clock; 
photons=[];  PE.n=0; photons.count=0; photons.sourceTime=0; photons.Gppid=0; photons.endID=0; photons.ancestor=0;
P=P(~all(cellfun('isempty', P),2),:);

p1                  = cat(1,P{:,1}); if isempty(p1); return; end
veluvecs            = cat(1,P{:,2});
t1                  = cat(1,P{:,3});
Gptid0              = cat(1,P{:,4});
wl                  = cat(1,P{:,5});

%SCINTILLATION PHOTONS IC
cherenkov0 = wl~=0;
i = ~cherenkov0; n=sum(i); %scintillation photons
if n>0
    veluvecs(i,:) = isovecs(n);
    t1(i) = t1(i) + randcdfc(input.Material(1).R.scintT,n);
    wl(i) = randcdfc(input.Material(1).R.scintWL,n);
end

if flags.status.fibers
    shapeID=2;  %shapeID = input.cube.shapeID; %[3-cylinder 2-cube]
    V = input.Volume(1);
    if shapeID==2 %cube
        ix = fcnindex1(V.vx,p1(:,1));  ix(ix<1 | ix>V.fn(1))=1;
        iy = fcnindex1(V.vy,p1(:,2));  iy(iy<1 | iy>V.fn(2))=1;
        v = sub2ind(V.fn,ix,iy);  v=min(v,size(V.xyz,1));
        dx = abs( V.xyz(v,[1 2]) - p1(:,[1 2]) );  r=sqrt(sum(dx.^2,2));
    else %arbitrary shape
        [v, r] = knnsearch(V.knnxyz,p1(:,[1 2]));
        dx = abs( V.knnxyz.X(v,:) - p1(:,[1 2]) );
    end
    uv = fcnunique(v);  n=numel(uv);  S=cell(n,1);
    
    if shapeID==3
        inside = r<V.L(1); %1mm fibers
    else
        inside = dx(:,1)<V.L(1) & dx(:,2)<V.L(2); %rectangular fibers
    end
    inside=find(inside); v=v(inside);
    
    [~,~,~,~,vx,vy,vz] = cube([V.L(1:2) input.cube.Lr(3)*1.001],[0 0 0],0); VX=cell(1,n); VY=VX; VZ=VX;
    for j=1:n
        x = V.xyz(uv(j),:);
        Lr = [V.L x];
        
        VX{j}=vx+x(1);  VY{j}=vy+x(2);  VZ{j}=vz;
        i = inside(v==uv(j)); %1mm fibers
        
        %N = [4 5 5 2]; %1.6-1.49 %fiber-cladding1
        %N = [4 6 6 2]; %1.6-1.42 %fiber-cladding2
        %N = [4 3 3 2]; %1.6-1.00 %fiber-air-pmt
        %N = [4 3 3 4]; %1.6-1.00 %fiber-air-fiber instead of pmt
        S{j} = transport1Volume(input.Material,[4 3 3 2],Lr,shapeID,p1,veluvecs,t1,wl,i);
    end
    if iscell(S{1}); S=cat(1,S{:}); end
    
    if n>0
        photons.fibervx = cell2mat(VX); %fiber rectangle
        photons.fibervy = cell2mat(VY);
        photons.fibervz = cell2mat(VZ);
    end
else
    S = transport1Volume(input.Material,[1 2 2 2],[input.cube.Lr 0 0 0],input.cube.shapeID,p1,veluvecs,t1,wl,(1:numel(i))'); 
end


try P=cat(1,S{:}); catch; P=S; end; if isempty(P); return; end
p1                  = cat(1,P{:,1});
t1                  = cat(1,P{:,2});
wl                  = cat(1,P{:,3});
endID               = cat(1,P{:,4});
t2                  = cat(1,P{:,5});
p2                  = cat(1,P{:,6});
ancestor            = cat(1,P{:,7});

Gptid = Gptid0(ancestor);
cherenkov = cherenkov0(ancestor);  %normal=1,  cherenkov=2


%FIND OUT WHICH PHOTONS PASS THE QE TEST ----------------------------------
af=input.cube.activefaces;  
if ~any(af==1 | af==3); endID(endID==11)=6; end
if ~any(af==2 | af==4); endID(endID==12)=6; end
if ~any(af==5 | af==6); endID(endID==13)=6; end


v3 = find(endID>10); %hit a wall and not reflected
nv3 = numel(v3);  nj=0;
if nv3>0
    QE = interp1c(input.wl, input.cube.QE, wl(v3));
    v4 = v3(rand(nv3,1)<QE);
    
    if numel(v4)>0    %ASSIGN PIXEL TO EACH PE
        ps = input.cube.pixelsize/2;
        pend = p2(v4,:);
        rpy = fcnVEC2RPY(-input.cube.all.normalVec);
        [j, r] = knnsearch(input.cube.all.knnxyz,pend);
        
        dv = rotateB2Wc(rpy(j,:), pend - input.cube.all.xyz(j,:));
        
        MTC766flag=false;
        if MTC766flag
            [~,QEmap]=fcnpruninglist;
            i = QEmap(j) & abs(dv(:,3))<ps(2) & abs(dv(:,2))<ps(1) & r<sqrt(sum(ps.*ps)) & abs(dv(:,1))<5;  k=v4(i);  j=j(i);  dv=dv(i,:);  %DELETE
        else
            i =            abs(dv(:,3))<ps(2) & abs(dv(:,2))<ps(1) & r<sqrt(sum(ps.*ps)) & abs(dv(:,1))<5;  k=v4(i);  j=j(i);  dv=dv(i,:);
        end
        
        %MICROCELL FILL FACTOR
        if flags.status.fibers && numel(j)>0
            fi=sub2ind(V.mcn,fcnindex1c(V.mcx,dv(:,2)),fcnindex1c(V.mcy,dv(:,3)));
            
            [X,si]=sortrows([j fi t2(k)],3);
            [~,ui]=unique(X(:,1:2),'rows'); sui=si(ui);
            j=j(sui);  k=k(sui);  dv=dv(sui,:);  %fprintf('Microcell Loss Rate: %.2f%%\n',100-numel(j)/numel(si)*100)
        end
        nj=numel(j);
    end
end

if nj>0
    PE.pixel        = j; %pixel ID
    PE.x            = p2(k,:); %pos in detector (mm)
    PE.t            = t2(k) + randcdfc(input.cube.pmt.transit.c,nj); %time (ns)
    PE.xpixel       = dv; %pos in pixel (mm)
    PE.amplitude    = randcdfc(input.cube.pmt.amplitude.c,nj);
    PE.Gptid        = Gptid(k);
    PE.Gpid         = fcntid2pid(PE.Gptid,G1.tid,G1.pid);
    endID(k)        = 7; %PE!
else
    PE.pixel        = [];
    PE.x            = [];
    PE.t            = [];
    PE.xpixel       = [];
    PE.amplitude    = [];
    PE.Gptid        = [];
    PE.Gpid         = [];
end

%DARK COUNTS
np = input.cube.pixels;
dArea = input.cube.pixelarea*np; %(mm^2) pixel surface area
dTime = input.cube.dsp.t(end); %(ns) listening time
nd = poissrndc(input.cube.pmt.darkCountRate*1E3*dArea*dTime*1E-9); %(80 kHz/mm^2) Dark Count Rate for SENSL J60035
if nd>0
    pid             = randi(np,[nd 1]);
    PE.pixel        = [PE.pixel;        randi(np,[nd 1])];
    PE.x            = [PE.x;            input.cube.all.xyz(pid,:)];
    PE.t            = [PE.t;            rand(nd,1)*dTime];
    PE.xpixel       = [PE.xpixel;       zeros(nd,3)];
    PE.amplitude    = [PE.amplitude;    randcdfc(input.cube.pmt.amplitude.c,nd).*darkCountCrossTalk(nd)];
    PE.Gptid        = [PE.Gptid;        zeros(nd,1)];
    PE.Gpid         = [PE.Gpid;         zeros(nd,1)];
end
PE.n = nj+nd;


%ASSIGN -------------------------------------------------------------------
photons.count               = numel(wl);
photons.sourcePos           = p1;
photons.sourceTime          = t1;
photons.Gppid               = fcntid2pid(Gptid,G1.tid,G1.pid);
photons.Gptid               = Gptid;
photons.cherenkov           = cherenkov;  %normal=1 or cherenkov=2
photons.wavelength          = wl;
photons.endID               = endID;
photons.endTime             = t2;
photons.endPos              = p2;
photons.ancestor            = ancestor;
photons.endIDlabel          = {'4. Re-emitted','5. Reflected','6. Transmitted','7. PE'};
%fprintf('Done (%.1fs).\n',etime(clock,startclock));
end


function c=darkCountCrossTalk(n)
%n=input number of dark count PEs
%c=output random cross-talk multiples
%x=[1 2 3 4]; y=[270 60 15 6]; F=fit(x', y'/sum(y), 'exp1');  %'poly1' or 'exp1' or 'rat01'
x=1:10; %y=F(x)./sum(F(x)); fig; stem(x,y,'filled')
cdfy=[0.774      0.94893      0.98846       0.9974      0.99942      0.99987      0.99998            1            1            1];
c = floor(fcnrandcdf(cdfy(:),x(:),n,'linear'));  %fig; histogram(c,'BinMethod','integers'); xlabel('SiPM Crosstalk Multiples')
end

function RQ = transport1Volume(Material,mi,Lr,shapeID,p1,veluvecs,t1,wl,i)
if isempty(i) 
    z3=zeros(0,3);  z=zeros(0,1);  RQ={z3,z,z,z,z,z3,z}; 
    return
end
RT =  Material(mi(1)).R.scintT;
RWL = Material(mi(1)).R.scintWL; RWL=RWL-RWL(1);
%PEsensor = [Material(mi).PEsensor];

p1                  = p1(i,:);
veluvecs            = veluvecs(i,:);
t1                  = t1(i);
wl                  = wl(i);    wl = max(min(wl,1998),83);
ancestor            = i;

T = [Material(mi(1)).X, Material(mi(2)).X(:,1), Material(mi(3)).X(:,1), Material(mi(4)).X(:,1)]; %X=[1-ir, 2-al, 3-at, 4-re, 5-scintillation]  SINGLE VOLUME
%[P1,T1,P2,T2,Status,Generation,WL,ancestorIndex]=photontransportcV3(p1',t1,veluvecs',wl,shapeID,Lr,T,RT,RWL);
[P1,T1,P2,T2,Status,Generation,WL,ancestorIndex]=photontransportcV4(p1',t1,veluvecs',wl,shapeID,Lr,T,RT,RWL);
RQ = {P1',T1,WL,Status,T2,P2',ancestor(ancestorIndex+1)};


% mg = 2000;  RQ=cell(mg,1);%max generations
% for g = 1:mg
%     np = numel(wl);
%     if np<3;  break;  end
% 
%     [p2,vel2,t2,status,wlnext,t1next]=photontransportc(p1,t1,veluvecs,int32(wl),shapeID,Lr,T,RT,RWL);
%     v5=find(status==4 | status==5); n5=numel(v5);    %5-reflected
%     
%     RQ{g} = {p1,t1,wl,status,t2,p2,ancestor};
% 
%     if n5>0 && g<mg %DEFINE RE-EMITTED PHOTON SOURCES
%         p1                     = p2(v5,:);
%         veluvecs               = vel2(v5,:); %really uvecs
%         t1                     = t1next(v5);
%         ancestor               = ancestor(v5);
%         wl                     = wlnext(v5);
%     else %no re-emissions
%         break
%     end
% end
% 
% 
% P=cat(1,RQ{:}); 
% p1                  = cat(1,P{:,1});
% t1                  = cat(1,P{:,2});
% wl                  = cat(1,P{:,3});
% endID               = cat(1,P{:,4});
% t2                  = cat(1,P{:,5});
% p2                  = cat(1,P{:,6});
% ancestor            = cat(1,P{:,7});
% 
% c=find(endID==13);  d=ancestor(c);  
% fig(1,3,1); dx=rangec(p2(c,:)-p1(d,:));  dt=t2(c)-t1(d);  v=dx./dt; c=v>0 & dt<15;  z={dx(c),dt(c),v(c)};
% sca; histogram(z{1},60); xlabel('PE dx (mm)');  
% sca; histogram(z{2},60); xlabel('PE dt (ns)');  
% sca; histogram(z{3},60); xlabel('PE v (mm/ns)');
% for i=1:3; sca(i); legend(sprintf('%.3g \\pm %.3g',mean(z{i}),std(z{i})),'Location','Best'); fcntight('y'); end


%FIBER PLOTS with V3
%c=find(Status==13);  r=rangec(P2'-P1');  
%R=accumarray(ancestorIndex+1,r); fig(1,1,'12cm'); histogram(R(R~=0),60); xlabel('PE dx (mm)'); axis tight;

end