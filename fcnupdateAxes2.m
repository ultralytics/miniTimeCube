function fcnupdateAxes2(input,handles,G1,photons)
h=handles.GUI.axes2; 
cla(h)
if photons.count==0; return; end

set(h,'NextPlot','add')
x = linspace(0,input.tMax.current,300);
v0 = photons.sourceTime<=input.tMax.current;

c = G1.upidcolor;
for i=find(G1.upidflag)
    times = photons.sourceTime(v0 & photons.Gppid==G1.upid(i));
    if ~isempty(times)
        n=hist(h,times,x); j=n~=0; n(j)=log(n(j)); %n=n/max(n);
        h2 = bar(h,x,n,1,'FaceColor',c{i},'EdgeColor',c{i});
        h3 = get(h2,'Children');
        set(h3,'FaceAlpha',.5)
    end
end

set(h,'XLim',[0 input.tMax.current]);
axis(h,'off','tight')
