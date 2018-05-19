function [p2,dt] = fcndetectorintercept(input,p1,vel)
%p1 is the vector origin, vel is the velocity vector
%endPos is the vector intercept on the detector wall, dt is the time it took to get there

%FIND OUT WHICH OF THE SIX PLANES IT HITS FIRST
if input.cube.shapeID==3 %cylinder
    if numel(p1)==1; repmat(p1,size(vel)); end

    r1 = input.cube.radius;
    p1xy = p1(:,[1 2]);
    velxy = vel(:,[1 2]);
    a = sum(velxy.^2,2);
    b = 2*sum(p1xy.*velxy,2);
    c = sum(p1xy.^2,2) - r1^2;
    dtHitRadius = quadratic(a,b,c);
    
    sz = sign(vel(:,3));
    dtHitPlane = (sz*input.cube.height/2 - p1(:,3))./ vel(:,3); %time needed to hit each of the 6 planes given initial pos and vel
    dt= min([dtHitRadius dtHitPlane],[],2); %time it takes each photon to strike the closest CCD plane

else
    vel(vel==0)=eps;
    sv = sign(vel);
    dtHitPlane = (sv.*input.cube.Lr - p1)./ vel; %time needed to hit each of the 6 planes given initial pos and vel
    dt = min(dtHitPlane,[],2); %time it takes each photon to strike the closest CCD plane
end
%[p2,dt,hitplane]=cubeinterceptc(p1,vel,input.cube.Lr);

p2 = p1 + vel.*dt; %position the photon hits a plane or attenuates, whichever comes first


