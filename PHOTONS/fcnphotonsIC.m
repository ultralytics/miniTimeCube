% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [P, G] = fcnphotonsIC(input, G, i)
ncs = 4; %number of cherenkov segments (between GEANT points)
% 1 photons.sourcePos
% 2 photons.sourceVel
% 3 photons.sourceTime
% 4 photons.Gppid
% 5 photons.Gptid
% 6 photons.cherenkov
% 7 photons.wavelength

i=find(G.tid==G.utid(i));  ni=numel(i);  P=cell((ncs+1)*ni,5);
outside1 = ~G.p1inside(i);
outside2 = ~G.p2inside(i);
tid = G.tid(i(1)); %trackid
t10 = G.t1(i);
t20 = G.t2(i);
de0 = G.de(i);
% pid = x(:,7); %particleid
% tid = x(:,8); %trackid
% ptid = x(:,9); %parent trackid
% G1.ke1 = x(:,18); %begin ke
% G1.ke2 = x(:,19); %end ke
% G1.de = x(:,20);
% G1.p1 = x(:,10:12);
% G1.p2 = x(:,13:15);
% G1.t1 = x(:,16);

if any(outside1)
    if all(outside1); return; end
    lii = find(outside1,1,'first');
    G.p1inside(i(lii:end)) = false;
    G.p2inside(i(lii:end)) = false;
    outside1 = ~G.p1inside(i);
    outside2 = ~G.p2inside(i);
    
    lit = t10(max(lii-1,1)); %last inside time
    ctid1 = fcntid2ctid(tid,G.tid,G.ptid); %children's track ids
    for j = 1:numel(ctid1)
        k = G.tid==ctid1(j);
        trackstarttime = min(G.t1(k));
        if trackstarttime>lit
            G.parentinside(k) = false;
        end
    end
end

if any(outside2)
    lii = find(outside2,1,'first');
    G.p2inside(i(lii:end)) = false;
    outside2 = ~G.p2inside(i);
end

if ~all(G.parentinside(i)) || all(de0==0)
    return
end

outside1or2 = outside1 | outside2;
pid = G.pid(i(1)); %particleid
ke10 = G.ke1(i);
ke20 = G.ke2(i);
p10 = G.p1(i,:);
p20 = G.p2(i,:);
vi = 1:ni;
Gparentinside = G.parentinside(i);

BQF=0;
flags.cherenkov=0;

j = find(G.upid==pid(1));
mass = G.upidmass{j};
charge = G.upidcharge{j};
name = G.upidname{j};

if mass>900; BQF=0.126; end
if charge~=0; flags.cherenkov=1; end
k=0; %k is P row index
for i = vi
    if ~Gparentinside(i); break; end
    de = de0(i);  if de==0;  continue;  end
    p1 = p10(i,:);
    p2 = p20(i,:);  dx = p2-p1;  r=sqrt(dx(1)^2+dx(2)^2+dx(3)^2);  %mm
    t1 = t10(i);
    t2 = t20(i);
    ke1 = ke10(i); %MeV
    ke2 = ke20(i); %MeV
    segmentoutside = false;
    
    if outside1or2(i) %if a particle leaves all it's later locations for it's track ID are no longer in the sim.
        segmentoutside=true;
    end
    
    devis = de;
    if BQF>0;  devis = fcnBQF(de,r,.126); end %de to visible de
    
    
%     %STRANGE DE/DX DIVERGENCES IN GEANT DATA
%     %[hf]=fig; hf.Tag='abc';
%     try
%         sca(evalin('base','handles.GUI.axes1'))
%         c=find(strcmp(name,{'proton','C9','e-','neutron'})); if isempty(0); c=5; end
%         %plot(findobj('Tag','abc'),de,devis,'.','Color',fcndefaultcolors(c))
%         
%         if strcmp(name,'proton') && r>.11
%             ''
%         end
%     catch
%         ''
%     end
    
    
    
    npe = input.Material(1).yield*devis; %number Poisson events
    if npe>1
        np = poissrndc(npe);
        if np>0
            k = k+1; %P index
            [spos, stime] = GEANTrandsegment(p1,p2,t1,t2,np);
            if segmentoutside
                vin = fcninsidedetector(spos, input);  np = sum(vin);  if np==0; continue; end
                spos = spos(vin,:);  stime=stime(vin,:);
            end
            
            zv = zeros(np,1);
            P{k,1} = spos;
            P{k,2} = zeros(np,3);
            P{k,3} = stime;
            P{k,4} = tid+zv;                    %photons.Gptid
            P{k,5} = zv;
        end
    end
    
    if flags.cherenkov && r>0 && fcnke2c(ke1,mass)>.5 %1/ir; %speed of light threshold to start cherenkov calculations (c)
        vke = linspace(ke1, ke2, ncs+1); %ncs = number of cherenkov segments (i.e. 4)
        speed = fcnke2c((vke(1:ncs)+vke(2:ncs+1))/2,mass); %c
        
        ct = input.Material(1).cherenkov;
        ic = fcnindex1c(ct.speed,speed(:));  npe = ct.npe(ic) *  (r/ncs);
        
        if sum(npe)<1; continue; end
        
        vf = linspace(0,1,ncs+1)';
        vt = linspace(t1, t2, ncs+1);
        vp = [p1(1)+vf*dx(1), p1(2)+vf*dx(2), p1(3)+vf*dx(3)];
        DCM_detector2particle = fcnVEC2DCM_B2W(dx)'; %dx = p2-p1
        for j = find(npe>1)
            %wlpdf = cherenkov_magnitude(input.wl, input.Material(1).X(:,1), speed(j)); %number of photons/mm
            %wlcdf = cumsum(wlpdf);
            %npe = wlcdf(end)*(r/ni);
            
            np = poissrndc(npe(j));
            if np>0
                k = k+1; %P index
                p1=vp(j,:); p2=vp(j+1,:);
                t1=vt(j);   t2=vt(j+1);
                [spos, stime] = GEANTrandsegment(p1,p2,t1,t2,np);
                if segmentoutside
                    vin = fcninsidedetector(spos, input);  np=sum(vin);  if np==0; continue; end
                    spos = spos(vin,:);  stime=stime(vin,:);
                end
                wl = fcnrandcdf(ct.wlcdf(:,ic(j)), input.wl, np, 'linear');
                ir = interp1c(input.wl, input.Material(1).X(:,1), wl);
                
                zv = zeros(np,1);
                P{k,1} = spos;
                P{k,3} = stime;
                P{k,4} = tid+zv;                    %photons.Gptid
                P{k,5} = wl;                        %photons.wavelength
                
                roll = rand(np,1)*(2*pi)-pi; %cherenkov photon roll angle
                
                %pitch = real( acos(1./(speed(j)*ir)) ); mean(pitch)*r2d
                %sinp=sin(pitch); %cherenkov photon pitch angle
                x = 1./(speed(j)*ir);
                sinp = sqrt(1 - x.^2);
                upf=zeros(np,3);  upf(:,1)=x;  upf(:,2)=sin(roll).*sinp;  upf(:,3)=cos(roll).*sinp; %uvecs particle frame
                P{k,2} = (upf*DCM_detector2particle) * speed(j); %uvecs in detector frame
            end %if np2>0
        end %for i=1:10
    end
    if segmentoutside; break; end
end

