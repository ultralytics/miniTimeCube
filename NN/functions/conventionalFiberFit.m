% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function xhat = conventionalFiberFit(I, T, threshold)
nw=size(I,2)/2;  vi{1}=1:nw; vi{2}=vi{1}+nw;
if nargin==2
    threshold = [0.006, 0.50];  % (mV, fraction)
end

%     %LAPPD
%     l=803; %length of anode, 200mm or 800mm;
%     speed = 190; %(mm/ns) pulse speed down anode
%     thr=[.006 .5];  method={'fraction','amplitude'};
%     for i=2
%         t=cell(1,4);  for j=1:4; [~,~,~,tj]=fcnpulsewidth(I(:,vi{j})',thr(i),[.006 1E3],[1 1E3],1E6,method{i});  t{j}=tj*.1; end
%         xa = speed*(t{1}-t{2})/2;  ta = (t{1}+t{2} - l/speed)/2; %(mm) and (ns)
%         xb = speed*(t{3}-t{4})/2;  tb = (t{3}+t{4} - l/speed)/2; %(mm) and (ns)
%         x=[xa xb ta tb];  z=T(:,[1 1 7 8]);  e=z-x; %std(e)    %fig(2,2,1.5); for i=1:4; sca; plot(x(:,i),z(:,i),'.'); end
%     end
%     speed=186;  l=1000;  xc=speed*(ta-tb)/2;  tc=(ta+tb - l/speed)/2; %(mm) and (ns)
%     x=[xa xb xc ta tb  tc];  z=T(:,[1 1 2 7 8 3]);  e=z-x;  xtof=x;  %fig(2,3,1.5); for i=1:6; sca; plot(x(:,i),z(:,i),'.'); end
%     std(e)

l=200; %length 200mm or 800mm, or 1000mm
speed = 78; %(mm/ns) group speed
method={'amplitude','fraction'};  xhat=cell(2,1);
for i=1:2 %1-FT 2-CFD
    t=cell(1,2);  for j=1:2; [~,~,~,t{j}]=fcnpulsewidth(I(:,vi{j})',threshold(i),[.0019 1E3],[10 1E3],1E6,method{i});  t{j}=t{j}*.2; end    
    xa = speed*(t{1}-t{2});  ta = (t{1}+t{2} - l/speed)/2; %(mm) and (ns)
    x=[xa ta zeros(size(ta))];  z=T(:,[1 2 4]);  e=z-x;  ci=t{1}~=0 & t{2}~=0;  x(~ci,:)=nan;
    %disp([mean(ci) std(e(ci,:))]); fig(1,3,1.25); for k=1:3; sca; plot(x(ci,k),z(ci,k),'.'); end; fcnmarkersize(1)
    xhat{i}=x;
end
