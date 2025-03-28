% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnRun(input)
clc; close all
% a=find(any(input.MTC.C{1,2},2)); %1089
% for i=1:size(input.MTC.C,1);  input.MTC.C{i,2}=input.MTC.C{i,2}(1122,:);  end

a=cell2mat(input.MTC.C(:,2));  b=a(:,4);  i=b~=0 & ~isnan(b);  b=b(i);  %save b.6.27.distanceAdjusted.mat b
%load b.6.27.mat;  
b=b(b>145 & b<180);  z=b;

fig(2,1,'19cm');
histogram(b,200,'DisplayName','Data');
[ya,xa]=histcounts(b,200); xa=xa(2:end)/2+xa(1:end-1)/2;  dxa=xa(2)-xa(1);

x0=[148.43         4.31       1.5789      0.99481]; %[mu sigma tf
%LB=[130 0 0 0]; UB=[155 20 20 10];
%[xhat, fx] = fmincon(@(x) fcnL(z,x), x0,[],[],[],[],LB,UB)
xhat=x0;

mu=xhat(1); tr=xhat(2); tf=xhat(3); s=xhat(4);
y=EES(xa,mu,tr,tf,s); y=y.*(mean(ya)/mean(y)); plot(xa,y,'linewidth',2,'DisplayName','Fit','Color',[1 1 1]*.7); xyzlabel('Time (bins)','Counts Per Bin'); axis tight; fcnfontsize(16)
sca; stem(xa,y-ya,'filled','DisplayName','Residuals'); xyzlabel('Time (bins)','Residuals'); fcntight('x equal');
xhat*.360

% y=EES(xa,mu,0.85/0.360,1.51/0.360,0.01/0.360);
% pulsewidth(y)*0.360.*dxa
end



function fx=fcnLEMG(z,x)
mu=x(1); s=x(2); tf=x(3);
x
y=fcnpdfEMG(z,mu,s,tf); %smeared
y=log(y+1E-323);
fx=-sum(y);
end


function fx=fcnL(z,x)
mu=x(1); tr=x(2); tf=x(3); s=x(4); x

y=EES(z,mu,tr,tf,s);

y=log(y);
fx=-sum(y);
end



function y=EES(z,mu,tr,tf,s)
minz=min(z); maxz=max(z); spanz=maxz-minz;

xa=linspace(minz-spanz*.1,maxz+spanz*.1,1000)';  dxa=xa(2)-xa(1);
y=exp(-(xa-mu)/tf)-exp(-(xa-mu)/tr);  y(xa<mu)=0;
ys = fcnpdfnorm(xa,xa',s)*y;  %smeared
ys=ys./sum(ys);
y=interp1c(xa,ys,z);
end