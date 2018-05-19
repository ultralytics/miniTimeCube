function []=fcnUpdate1ParticleType(input, upidhandle, G, i)
pid = G.upid(i); %particleid
v0 = G.t1<input.plotTime;
v1 = find(G.pid==pid & v0);
if isempty(v1)
   upidhandle.XData=[];  upidhandle.YData=[];  upidhandle.ZData=[];
   return
end

utid = fcnunique(G.tid(v1));  nutid = numel(utid);
p3x=[]; p3y=[]; p3z=[];  nan13=nan(1,3);
for j=1:nutid
    tid = utid(j); %trackid
    v1=find(G.tid==tid & v0);
    
    %ptid = G.ptid(v1(1)); %parent trackid
    %ke1 = G.ke1(v1); %parent trackid
    %ke2 = G.ke2(v1); %parent trackid
    %de = G.de(v1); %parent trackid
    p1 = G.p1(v1,:);
    p2 = G.p2(v1,:);
    t1 = G.t1(v1);
    t2 = G.t2(v1);
    
    if t2(end)>input.plotTime
        p2=qinterp1([t1(end) t2(end)],[p1(end,:); p2(end,:)],input.plotTime);
    else
        p2=p2(end,:);
    end
    
    
    if nutid>1
        p = [p1; p2; nan13];
        
        p3x = [p3x; p(:,1)];
        p3y = [p3y; p(:,2)];
        p3z = [p3z; p(:,3)];
    else
        p = [p1; p2];
        
        p3x = p(:,1);
        p3y = p(:,2);
        p3z = p(:,3);
    end

end
%h=plot3(p3x,p3y,p3z,'.-','color',c,'MarkerSize',15);
upidhandle.XData=p3x;  upidhandle.YData=p3y;  upidhandle.ZData=p3z;

    
    
