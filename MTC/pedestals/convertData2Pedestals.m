% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = convertData2Pedestals(input)
A=input.MTC.A; input.MTC.A=[]; fname=A.filename;
E = A.E-min(A.E)+1;
W = A.x(:,8)+1;
pid = MTCSRCCH2pixelID(A.x(:,3:6));

ne=max(E);  X0=zeros(64,512,ne,'single');  n=cell(ne,1);  mu=n;  sigma=n;  Ax=A.x(:,10:73); A.x=[];
ci = 1:1536;
waitbarParfor(numel(ci));
for i=ci
    fprintf('%s\n',waitbarParfor);  X=X0;
    for l=find(pid==i)'
        X(:,W(l),E(l)) = Ax(l,:);
    end
    X=reshape(X,[64*512 ne]);
    [sigma{i}, mu{i}, n{i}] = nonzerostd(single(X),2);
end
mu=[mu{:}]';  sigma=[sigma{:}]';  n=uint16([n{:}]');
save(sprintf('PedestalData.%s.mat',fname),'mu','sigma','n','fname','ci');


%FIND OUTLIERS
muvb = int16(fcnchannels2windows(mu)'); %#ok<*NASGU>

[~,ia]=fcnsigmarejection(double(mu),6,3);
[~,ib]=fcnsigmarejection(double(sigma),6,3);  
inliers=ia & ib & mu>0 & sigma>0 & n>2;

%SAVE
mmap=load('MTCmapping.mat');
pruneChannels=mean(inliers,2)<.9;    badFraction=1-mean(inliers,2); badFraction=badFraction(mmap.X(:,1)); %FOR MAP
pedestalOutliers = fcnchannels2windows(~inliers)';  pedestalOutliers(:,sum(pedestalOutliers)>32)=true;
save('-v6',sprintf('pedestals.%s.mat',fname),'muvb','pedestalOutliers','pruneChannels');


fig(3,2,1.5); if numel(ci)==1536; cj=mmap.X(ci,1); else cj=1:numel(ci); end
cl={[1800 3000],[0 16],[0 30]};  cm={'redblue','bluewhite','bluewhite'};  x={mu,sigma,n};  s={'Mean','Sigma','Samples'};
for i=1:3
    b=x{i};  b(b==0)=nan;
    sca;    a=b(cj,:);       imagesc(a(:,1:1:end),cl{i});  colormap(gca,cm{i});  title(s{i});    axis image normal xy off; colorbar;
    h=sca;  a=nanmean(b,2);     a(a==0)=nan; fcnplotdetectorprojection(input,ci,a+.02,[0 0 0],'winkeltripel','SCROD'); colormap(gca,cm{i}); h.CLim=cl{i};  if i==1; title(str_(fname)); h.Title.Visible='on'; end; h.CameraViewAngle=h.CameraViewAngle*.85;
    %sca;    %a=double(b(inliers));    li=cl{i}; xh=linspace(li(1),li(2),100);  fhistogram(double(b),xh);     fhistogram(a,xh);   xyzlabel(s{i}); legend show
    %if i==3; h=gca; h.YScale='log'; end
end
%fcntight

fig(2,1,1,2); 
sca; imshow(inliers(mmap.X(:,1),:),[0 1],'XData',[1 512]); axis image normal xy; colormap(gray); title(sprintf('%.2f%% Outliers\n%s',(1-mean3(inliers))*100,str_(fname))); fcnfontsize(10)
h=sca; a=mean(inliers,2); fcnplotdetectorprojection(input,ci,a-.02,[0 0 0],'winkeltripel','SCROD'); colormap(gca,gray); h.CLim=[0 1]; h.CameraViewAngle=h.CameraViewAngle*.8;

''



