% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = fcnNeutronExits(input,output,flags,PE,G1,handles,plotFlag)
%x=[G1.p1inside G1.p2inside G1.parentinside];
i=find(G1.p1inside & G1.parentinside & G1.tid==1,1,'last');
results.xhat=zeros(1,29); results.true=zeros(1,29);

exit.E=G1.ke1(i);
if G1.p2inside(i)
    exit.pos=G1.p2(i,:);
else
    vec=G1.p2(i,:)-G1.p1(i,:);
    exit.pos = fcndetectorintercept(input,G1.p1(i,:),vec);
end

results.true(1:6)=[exit.pos exit.E G1.ke1(1) G1.p2inside(i)];
end


function plotTSresults
T=MC.xtrue;
i=T(:,1)==25;
x=T(:,5); %E0 (MeV)
fig; sca(1); [mu,sigma,~,p] = movingMean(x,i,100,0);            plot(p,mu,'-','Display',sprintf('\\mu=%.3g',mean(sigma)));  xyzlabel('Neutron E_0 (MeV)','Exiting Opposite (fraction)')

fig; plot3(T(:,1),T(:,2),T(:,3),'.'); axis equal; fcnmarkersize(1); fcnview('skew')

end