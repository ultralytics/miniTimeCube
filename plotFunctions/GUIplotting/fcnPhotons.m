function x = fcnPhotons(input,handles,flags,photons,PE,G,i)
v1 = find(photons.started & photons.Gppid==G.upid(i));

%GUI TABLE ----------------------------------------------------------------
x = cell(1,6);
%ended = photons.endTime < input.plotTime;
%started = photons.sourceTime < input.plotTime;
%alive = started & ~ended;
j = find(G.pid==G.upid(i));  nPID=numel(j);
x{2} = G.nupidv(i); %nTID
x{3} = nPID;
x{4} = double(sum(G.de(j)));
x{5} = numel(v1);
if PE.n>0; x{6} = sum(PE.t(PE.Gpid==G.upid(i))<input.plotTime); else x{6}=0; end


if isempty(v1)
    return
end


%PLOT PHOTONS -------------------------------------------------------------
size = max(input.cube.Lr)*.02 * input.photonsize/4;

if flags.status.plotphotonsalive
    v4 = v1( photons.endTime(v1)>input.plotTime); %still alive
    if numel(v4)>0  %number alive
        t1 = photons.sourceTime(v4,1);

        a=photons.sourcePos(v4,:); b = photons.endPos(v4,:);  
        [r,dx]=rangec(a,b);
        
        f = (input.plotTime - t1)./(photons.endTime(v4,:) - t1);  pos = a + f.*dx;
        
        plotPhotons(dx,f,r,pos,v4,size,photons,handles,i);
    end
end

if flags.status.plotphotonsource
    pos = photons.sourcePos(v1,:);
    plotPhotons(pos,v1,size,photons);
end

if flags.status.plotphotonend
    v6 = v1(photons.endTime(v1)<input.plotTime & photons.endID(v1)==6); %transmitted!
    
    if numel(v6)>0
        pos = photons.endPos(v6,:);
        plotPhotons(pos,v6,size,photons);
    end
end


end

function plotPhotons(dx,f,r,pos,i,size,photons,handles,hi)
dxu = dx./r;
p1 = pos - dxu.*min(f.*r,size);

A = nan(3,numel(i));
j = photons.cherenkov(i);  k=~j;
hc=handles.upidcp(hi);
hs=handles.upidsp(hi);
A(1,:) = p1(:,1); A(2,:)=pos(:,1);  a=A(:,k); hs.XData=a(:);  a=A(:,j); hc.XData=a(:);
A(1,:) = p1(:,2); A(2,:)=pos(:,2);  a=A(:,k); hs.YData=a(:);  a=A(:,j); hc.YData=a(:);
A(1,:) = p1(:,3); A(2,:)=pos(:,3);  a=A(:,k); hs.ZData=a(:);  a=A(:,j); hc.ZData=a(:);
end

