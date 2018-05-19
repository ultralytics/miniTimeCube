function xa=prepareRCDF(x,c) %(x,cdf)
%this function prepares a cdf for random value lookup by randcdfc.cpp
%example:  fig; fcnhist(randcdfc(xa,5E4),100);
x=x(:);
c=c(:);
ca=linspace(0,1,16384);


i = c>circshift(c,1);
if sum(i)==1; j=find(i,1,'first'); i(j-1)=true; end %if only 1 unique value, i.e. cdf=[0 0 1 1 1...]
if any(i)
    d=c(i);  d=d-min(c);  d=d/max(d);
    F = griddedInterpolant(d,x(i),'linear');
    xa=F(ca);
else
    xa=zeros(1,16384);
end



