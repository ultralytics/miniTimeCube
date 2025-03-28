% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [input,output,G1,photons,PE,flags] = fcn1Event(input,flags,handles,G,V) 
photons.count = 0;
G1 = [];
timeOffset = 25;%rand*3+45; %(ns), TYPICALLY ZERO
energyGain = 1; %TYPICALLY 1

if exist('V','var') && ~isempty(V)
    input.eventNumber = V(4);
    timeOffset = V(2);
    energyGain = V(3);
else
    input.eventNumber = eval(handles.GUI.edit1.String);
end
C_G2N = [   0 0 1
            0 1 0
            1 0 0   ]; %GEANT 2 NED
%DEFINE PIXEL INFO --------------------------------------------------------
[input, flags] = fcnDefinePixels(input,flags,handles);
flags.update.detectorgeometry=0; %geometry all updated

%INCOMING ANGLE ----------------------------------------------------
if flags.status.vertexfixedangle %if handles.GUI.checkbox59.Value==1
    vel = fcnvec2uvec([-1 0 0]); %fixed
else
    vel = isovecs(1); %random
    %vel = diffusevecs(1,[0 0 1]);
end
rpy = fcnVEC2RPY(vel);  %rpy(1) = pi*(rand*2-1);
C = C_G2N*fcnRPY2DCM_W2B(rpy);

%RANDOM INSIDE VOLUME
if flags.status.vertexcentered
    offset = [2.5 2.5 0]; %centered
else
    switch input.cube.shapeID
        case 3 %cylinder
            cc = fcnSC2CC(input.cube.radius*sqrt(rand),  0,  rand*2*pi-pi);
            offset = [cc(1:2) input.cube.height*(rand-1/2)];
        otherwise %cube
            offset = 2*input.cube.Lr.*(rand(1,3)-.5);
            if flags.status.fibers && flags.status.MC && strcmp(G.particleName,'PointSource')
                offset(1:2) = 2.5+(rand(1,2)-.5)*5; %5mm FIBERTRAIN 64CH
                %offset(3) = -421;
            end
    end
end

%FROM A FIXED POINT
% particleOrigin = [5000 0 0];
% vel = fcnvec2uvec([offset - particleOrigin]);
% rpy = fcnVEC2RPY(vel);  rpy(1) = pi*(rand*2-1);
% C = C_G2N*fcnRPY2DCM_W2B(rpy);


%RANDOM ON DETECTOR WALL
if flags.status.vertexonwall %if handles.GUI.checkbox52.Value==1
    p1 = offset;  vel = -vel;
    if input.cube.shapeID==3 %cylinder
        p1xy = p1(:,[1 2]);
        velxy = vel(:,[1 2]);
        a = sum(velxy.^2,2);
        b = 2*sum(p1xy.*velxy,2);
        c = sum(p1xy.^2,2) - input.cube.radius^2;
        dtHitRadius = (-b+sqrt(b.^2-4*a.*c))./(2*a); %quadratic
        
        p1z = p1(:,3);
        velz = vel(:,3);
        dtHitPlane = (sign(velz)*input.cube.height/2 - p1z)./ velz; %time needed to hit each of the 6 planes given initial pos and vel
        dt = min([dtHitRadius dtHitPlane],[],2); %time it takes each photon to strike the closest CCD plane
    else
        dtHitPlane = abs((sign(vel).*input.cube.Lr - p1)./ vel); %time needed to hit each of the 6 planes given initial pos and vel
        dt = min(dtHitPlane,[],2); %time it takes each photon to strike the closest CCD plane
    end
    offset = p1 + vel*dt;
end
G1.p0 = offset;


%GEANT STUFF --------------------------------------------------------------
en = mod(input.eventNumber,G.ne); %eventnumber
i = G.ve==en;
v1 = G.ei(i,1):G.ei(i,2);
if ~any(v1) %NO GEANT DATA
    [input, output] = fcnBuffer2Memory(input,flags,handles,photons);
    G1.na = 0;  G1.nua = 0;  G1.uake1 = 0;
    G1.nupid = 0;  G1.upid = [];
    G1.nutid = 0;  G1.utid = [];
    G1.uvec(1,:) = vel./norm(vel);
    G1.upidflag = ones(1,G1.nupid);
    G1.nupidv = zeros(1,G1.nupid);
    input.tMin = 0;
    input.tMax.current = 100;
    fprintf('Warning!! No GEANT data found for this event.\n')
    return
end
x = double(G.x(v1,:));
if isfield(G,'double')
    x(:,2) = G.double(v1); %PDGEncoding has to be read as a double precision to not lose accuracy
end

switch size(x,2)
    case 10 %NEW MITCH 2015
        x = sortrows(x,[3 8]); %sort by tid then time
        pos=x(:,5:7)*C + offset;
        j=[false; diff(x(:,3))==0]; %find changes in tid
        if x(1,2)==-12 %nuebar
            j(1)=true;
            a=max(x(:,3))+1; x(1,3)=a; %new nuebar trackID
            x(x(:,4)==0,4)=a; %set n and e+ to have this parent id
            x(1,4)=0; %nuebar parent id = 0
            j=find(j); %post-step
            i=j-1; i=max(i,1); %prestep
        else
            j=find(j); %post-step
            i=j-1; %prestep
        end
        
        %a=x(j,10)>0; if any(a); t0=min(x(i(a),8)); x(:,8)=x(:,8)-t0; end %OPTIONAL. Sets first dE time to 0 ns
        
        G1.p1  = pos(i,:);
        G1.p2  = pos(j,:);
        G1.t1  = x(i,8)+timeOffset;
        G1.t2  = x(j,8)+timeOffset;
        G1.ke1 = x(i,9)*energyGain;
        G1.ke2 = x(j,9)*energyGain;
        G1.de  = x(j,10)*energyGain;
        
        pid    = x(j,2); %particleid
        tid    = x(j,3); %trackid
        ptid   = x(j,4); %parent trackid
        nx=numel(j);
    otherwise %OLD
        x = sortrows(x,[8 17]); %sort by tid then time
        G1.p1  = x(:,10:12)*C + offset;   nx=size(x,1);
        G1.p2  = x(:,13:15)*C + offset;
        G1.t1  = x(:,16)+timeOffset;
        G1.t2  = x(:,17)+timeOffset;
        G1.ke1 = x(:,18)*energyGain; %begin ke
        G1.ke2 = x(:,19)*energyGain; %end ke
        G1.de  = x(:,20)*energyGain;

        pid    = x(:,7); %particleid
        tid    = x(:,8); %trackid
        ptid   = x(:,9); %parent trackid
end
G1.pname = fcnpid2name(pid);
G1.p1inside = fcninsidedetector(G1.p1, input);  G1.p1inside(1)=true;
G1.p2inside = fcninsidedetector(G1.p2, input);
G1.parentinside = true(nx,1);
G1.x = x;


G1.pid = pid; %particle ID
G1.tid = tid; %track ID
G1.ptid = ptid; %parent track ID
G1.upid = fcnunique(pid)';  G1.nupid=numel(G1.upid);
[G1.utid, G1.utidi] = fcnunique(tid);  G1.nutid=numel(G1.utid); %utidfi = unique track ID's first indices
G1.ppid = fcntid2pid(ptid,tid,pid); %parents pardicle ID


%DENDROGRAM ANCESTRY MATRIX -----------------------------------------------
x = nan(G1.nutid,30);  x(:,1) = G1.utid;
J = zeros(max(G1.utid),1); J(G1.utid) = 1:G1.nutid;
for i = 2:30 %generations
    j = x(:,i-1)>0;
    if ~any(j); break; end
    x(j,i) = G1.ptid(G1.utidi(J(x(j,i-1))));
end
x = x(:,1:i-1);


G1.dendrogram = x;
minx = minnonzero(x,2);
G1.atid = minx(J(G1.tid)); %ancestor track ID
G1.ua = fcnunique(G1.atid); %unique ancestors (ptid==0)
G1.nua = numel(G1.ua); %number of ancestors

%ANCESTOR UVECS -----------------------------------------------------------
G1.uake1 = zeros(G1.nua,1);
G1.uvec = zeros(G1.nua+1,3);
G1.uvec(1,:) = [0 0 1]*C; %always the same angle
for i=1:G1.nua
    j=find(tid==G1.ua(i),1,'first');
    G1.uvec(i+1,:) = fcnvec2uvec(G1.p2(j,:)-G1.p1(j,:));
    G1.uake1(i) = G1.ke1(j); %unique ancestor initial ke1
end
if norm(G1.uvec(1,:))==0; G1.uvec(1,:)=G1.uvec(2,:); end

G1.upidflag = ones(1,G1.nupid);
G1.upidv{i} = cell(1,G1.nupid);
G1.nupidv = zeros(1,G1.nupid); 
for i = 1:G1.nupid
    v = fcnunique(G1.tid(G1.pid==G1.upid(i)));
    G1.upidv{i} = v; %utids per upid
    G1.nupidv(i) = numel(v); %number of utids per upid
end
[G1.upidname,G1.upidcolor,G1.upidmass,G1.upidcharge] = fcnpid2name(G1.upid);


%GET PHOTONS --------------------------------------------------------------
P = cell(G1.nutid,1);
for i=1:G1.nutid
    [P{i}, G1] = fcnphotonsIC(input, G1, i);  
end
P = cat(1,P{:});
[photons, PE]= fcnphotontransport(input, P, G1, flags);



%FIND CURRENT TMAX --------------------------------------------------------
input.tMin = 0;
if photons.count>0;  input.tMax.current=max(photons.endTime)+10;  else  input.tMax.current=max(G1.t2);  end


% if flags.status.MC
%    n=V(1);
%     %n=50;
%     i=find(PE.t<25 & PE.pixel==1 & abs(PE.xpixel(:,3))<90); ni=numel(i);
%     PE.pixel(:)=0; PE.pixel(i(randperm(ni,n)))=1;
% end

%DUMP MCP BUFFER MEASUREMENTS INTO MEMORY ---------------------------------
[input, output] = fcnBuffer2Memory(input,flags,handles,PE);

if photons.count>500000
    fprintf('WARNING! LARGE PHOTON COUNT: %.2fM photons created\n',double(photons.count)/1E6)
end

