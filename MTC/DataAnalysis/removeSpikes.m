function V=removeSpikes(V,threshold)
if nargin==1; threshold=50; end

if numel(threshold)==1 %REMOVE SINGLE SPIKES
    a1 = circshift(V(:),-1);
    a3 = circshift(V(:), 1);
    
    r13 = abs(a1-a3)*3;
    r12 = abs(a1-V(:));  clear a1
    r32 = abs(a3-V(:));  clear a3
    
    j = (r12>threshold & r32>threshold) & (r13<r12 | r13<r32);  clear r12 r13 r32
    
    j = find(j);  j=j(rem(j,64)~=0);   n=numel(V);
    if any(j<2 | j>(n-1)); j=floorandceil(j,2,n-1); end
    if numel(j)>0;  V(j)=(V(j-1)+V(j+1))/2;  end
    
else %numel(threshold)==2; %REMOVE STATISTICAL OUTLIERS
    ec = ~any(diff(V)); 
    if any(ec) %if any empty channels, do not operate on them
        if all(ec); return; end
        V0=V; V=V(:,~ec); 
    end 

    srl=threshold(1); ni=threshold(2);  %srl=3; ni=3;
    dx=fcndiff(V(:),1,1);
    [~,ix]=fcnsigmarejection(dx,srl,ni,'onlyIndices');  clear dx
    dv=fcndiff(V(:),2,1);  
    [~,iv]=fcnsigmarejection(dv,srl,ni,'onlyIndices');  iv=circshift(iv,-1); clear dv
    ys = ~(ix & iv);  if ~any(ys); if any(ec); V=V0; end; return; end  %ys = outliers
%    xv=1:1024;  fig; plot(xv,V,'.-'); plot(xv(~ix),V(~ix),'ro'); plot(xv(~iv),V(~iv),'ro'); fcntight
    yt=(circshift(ys,1) | circshift(ys,-1)) & ~ys;     %yt=fcnsmooth(ys,3) & ~ys; %good points
    i=find(yt);  F=griddedInterpolant(i,V(i),'linear','none'); %DOUBLES NEEDED HERE!
    j=find(ys);  V(j)=F(j);

    if any(ec); V0(:,~ec)=V; V=V0; end
end





