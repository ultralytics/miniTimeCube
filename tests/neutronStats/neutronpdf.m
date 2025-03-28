% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function neutronpdf(input,flags,handles,GEANT)
clc; tic
r2d = 180/pi;
nj=10000; %GEANT.ne;
pname = fcnpathm1(mfilename('fullpath'));
%x = [1-pid, 2-neutron range (mm), 3-time (ns), 4-angle (rad), 5-dE (MeV), 6-emitted photons, 7-observed photons, 8-inside detector, 9-x, 10-y, 11-z, ...
%     12-xhat, 13-yhat, 14-zhat, 15-dEhat] %PER BOUNCE
%xcapture = [neutroncapturexyz, 4-capturetime, 5-capturephotons, 6-allinsideflag, 7-E0, 8-E0enclosed, 9-E0hat]
%2212 = proton
%1000060120 = C12

ni = numel(GEANT.files.names);
flags.status.MC = 1;
x = zeros(nj,15,2,ni); %neutron bounce info
xcapture = zeros(nj,9,ni); %neutron capture info
for i=1:ni
    GEANT = fcnloadGEANT(GEANT.files.names{i});  
    for j=1:nj
        fprintf('%.0f/%.0f\n',j,nj)
        [input1,output1,G1,~,PE] = fcn1Event(input,flags,handles,GEANT,[0 0 0, j-1]);
        
        E0hat = 0;
        %if output1.zN.sum(1)>2
        %    xout = fcnMLpoint(input,output1.zN.z{1}); E0hat = xout(4);
        %end
        k = find(G1.tid==1,1,'last');
        E0 = G1.ke1(find(G1.tid==1,1,'first'));
        E0enclosed = sum(G1.de(G1.t1<35 & G1.parentinside & G1.p1inside & G1.p2inside & G1.pid~=1000030070 & G1.pid~=1000020040));
        xcapture(j,:,i) = [G1.p2(k,:) G1.t2(k) output1.zN.sum(2) G1.parentinside(k) E0 E0enclosed E0hat]; %#ok<*PFOUS>
        
        nv = find(G1.tid==1); %neutron rows
        for b=1:2 %bounce
            k = find(G1.tid==(b+1),1,'first'); %2 = first bounce, 3 = second bounce, etc.
            if numel(k)>0 && numel(nv)>b
                
                p2 = G1.p2(nv(b),:);
                vec1 = G1.p2(nv(b),:)   - G1.p1(nv(b),:); %neutron vec before bounce
                vec2 = G1.p2(nv(b+1),:) - G1.p1(nv(b+1),:); %neutron vec after bounce
                
                v3 = PE.Gptid==(b+1);
                range = norm(vec1); %mm
                time = G1.t2(nv(b)) - G1.t1(nv(b)); %ns
                angle = fcnangle(vec2,vec1); %rad
                dE = G1.ke1(nv(b)) - G1.ke2(nv(b)); %MeV
                emittedPhotons = sum(v2);
                observedPhotons = sum(v3);
                insidedetector = fcninsidedetector(p2, input1);
                
                %RUN CHOOZ ESTIMATOR
                xout = [0 0 0 0];
                %if sum(v3)>2
                %    zN = accumarray(photons.endPixelID(v3), 1, [input1.cube.pixels 1]);
                %    xout = fcnMLpoint(input,zN);
                %end
                x(j,:,b,i) = [G1.pid(k), range, time, angle*r2d, dE, emittedPhotons, observedPhotons, insidedetector, p2, xout(1:4)];
            end
        end
    end
end
toc
save(fcnincrementfname([pname 'NewGEANT 135MeV Neutron Statistics.mat']),'x','xcapture')

%i = 5; neutronpdfplots(input,G1,handles,x(:,:,:,i),xcapture(:,:,i))

