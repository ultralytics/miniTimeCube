function V = subtractpedestals(pedestals,pedestalOutliers,pid,window,V)
% load('MTCpedestal') %loads X 98304x71
% pedestals = muvb(:,8:end);  %#ok<NODEF> %idxm=load('MTCmapping.mat');  b=b(idxm.X(:,1),:);

n=numel(V);
% if mean(V(1:min(n,1E4)))<1000 %pedestals already subtracted
%     V(V<-1200)=1600; 
%     return
% end


nw=512; %64 (number of windows)

j=V<100;
i=sub2ind([nw 1536],window+1,pid);
V = V - pedestals(:,i);

%MTC ONLY!!
V(j)=1600;
%V(pedestalOutliers(:,i))=1600; %not nan because int32

