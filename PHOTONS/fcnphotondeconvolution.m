% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [V, QT] = fcnphotondeconvolution(input,flags,X)
%signal = 1536x256
if isempty(X); V=X; QT=[]; return; end
firstphotonsonly = true;
pmt = input.cube.pmt;  dsp = input.cube.dsp; 
sX = size(X);
X = reshape(X,[sX(1) prod(sX(2:end))])';

%PSF
p=[];
p.t = dsp.dt*128;
p.pixel = 1;
p.xpixel = [0 0 0];
psf = squeeze(fcnanalogvoltage(input,flags,p,dsp.t,0,1,false))';  edge = 0;  psf = psf(1+edge:numel(psf)-edge);
psf = psf(:)'/sum(psf);


%WIENER DECONVOLUTION
nsr = double( input.cube.dsp.noise^2/var(X(:)) ); %4mV/var(X)
if flags.status.lappd    
    V = deconvwnr(X,psf,nsr);  %fig; plot(V')
    threshold=.1;
else
    V = deconvwnr(X,psf,nsr*1);
    threshold = .1; %VERY IMPORTANT!!!!
end
[pixeli, timei, maxval] = findlocalmaximaregions(V',threshold); if ~any(pixeli); QT=[]; return; end


%pulse widths (bins)
Xt = X';
width = zeros(size(Xt,2),1);

corrections = width*0;
QT = [pixeli,  timei(:)-corrections(pixeli),  width(pixeli)];


if firstphotonsonly
    [~,j] = fcnunique(pixeli); %photons already sorted by index first then time
    QT = QT(j,:);
end
QT = QT(all(isfinite(QT),2),:);

if flags.status.lappd;  [stripi,sidei,tilei] = ind2sub(sX(2:end), pixeli);  end %#ok<ASGLU>

end

function [cols, tmax, maxval] = findlocalmaximaregions(y,threshold)
%y = nxm, finds the maxima down dimension x above threshold
cols=[]; tmax=[]; maxval=[]; sy=size(y);
if nargin<2;  threshold=150;  end;  y0=y;  y(y<threshold)=0;

dy = fcndiff(y,1);  dsdy = fcndiff(sign(dy),1);
ind = find(dsdy==-2) - 1; %falling zero crossings

row = mod(ind,sy(1));  ind=ind(row>2 & row<(sy(1)-2));  ni=numel(ind);
if ni==0; return; end
    
[i, cols] = ind2sub(size(y),ind); %i=row, j=col
i = bsxfun(@plus,i,-2:2);
ind = bsxfun(@plus,ind,-2:2);

%UPSAMPLE
nss = 300; %number spline samples
xi = linspace(1,3,nss);
yi = interp1((0:4)',y0(ind)',xi,'spline');  if ni>1; yi=yi'; end
xip = bsxfun(@plus, xi, i(:,1));

[maxval, j]=max(yi,[],2);
tmax = xip(sub2ind(size(yi),1:numel(j),j'));

%PLOT
%y(y<threshold) = nan;  x=(1:size(y,1))';
%fig; plot(x,y,'.-'); plot(i,y(ind),'b.','markersize',20); plot(xip,yi,'r.','Markersize',5); plot(tmax,maxval,'mo')
end
