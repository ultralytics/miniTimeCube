% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function fcndeconvolve(input,output,flags)
clc; close; 
pmt = input.cube.pmt;  dsp = input.cube.dsp;
h=fig(2,2,'19x38cm'); 

%SIGNAL -------------------------------------------------------------------
ni = 3; %number of PE to make
nt = dsp.samples;
PE = [];
%PE.t = rand(ni,1)*45;
PE.t = randcdfc(input.Material(1).R.scintT,ni)+20;
PE.pixel = ones(ni,1);
PE.amplitude = randcdfc(input.cube.pmt.amplitude.c,ni);
V = fcnanalogvoltage(input,flags,PE,dsp.t,0,1:ni,true);

pid = PE.pixel(1:ni);
[~, I, J] = unique(pid,'rows');  nus = numel(I);  pid = pid(I);
c = fcnaccumrows(V,J);

n = 1;
c = permute(c,[2 1]);
X = zeros(nt, n);
X(:,pid) = c;


MTCdataflag=true;
if MTCdataflag
    D = output(1).D;
    i=any(D>0,1);
    D = D(:,i);  
    %=randi(size(D,2),1);
    [~,j]=sort(max(D,[],1));
    signal = D(1:256,j(end-55))';
    
    %PSF ----------------------------------------------------------------------
    i=1;
    p2=PE;
    p2.t(i) = dsp.dt*128;
    p2.pixel(i) = 1;
    psf = fcnanalogvoltage(input,flags,p2,dsp.t,0,i,false);  edge = 0;  psf = psf(1+edge:256-edge);  psft = dsp.t(1+edge:256-edge);
    nsr = 13^2 / var(signal);
    dsp.t = 1:size(signal,2);
    threshold=50;
else
    %NOISE --------------------------------------------------------------------
    B = 25*dsp.dt; %autocorrelation time constant
    noise = randncorr([n nt], B/dsp.dt)*dsp.noise;
    X = X + reshape(noise',size(X));
    X = floorandceil(X,dsp.dynamicrange);
    
    nb = 2^dsp.bits;
    X = round(X*nb)/nb;
    signal = X';
    
    V2 = V; V2(V2==0)=nan;
    sca(h(1)); plot(dsp.t,V2,'.-'); title(sprintf('Signal (%.0f PE)',ni));
    sca(h(3)); plot(dsp.t,noise,'.-'); title(sprintf('Noise (\\tau = %.2fns, \\mu = 0mV, \\sigma = %.1fmV)',B,dsp.noise*1e3));
    
    %PSF ----------------------------------------------------------------------
    i=1;
    p2=PE;
    p2.t(i) = dsp.dt*128;
    p2.pixel(i) = 1;
    psf = fcnanalogvoltage(input,flags,p2,dsp.t,0,i,false);  edge = 0;  psf = psf(1+edge:256-edge);  psft = dsp.t(1+edge:256-edge);
    nsr = input.cube.dsp.noise^2 / var(signal);
    threshold = max(signal)*.2;
end
psf = psf/sum(psf); %MUST SUM TO UNITY
sca(h(2)); 
plot(dsp.t,signal,'.-','Display','Signal + Noise'); 
%plot(psft,psf,'g.-','Display','PSF');
fcntight(h(1:4),'jointy'); legend show


%WIENER DECONVOLUTION -----------------------------------------------------
sca(h(4)); tic

%icorr = autocorr(psf(55:end),52);
%ncorr = exp(-dsp.t/B)/1E4;
%y = deconvwnr(pdf,psf,ncorr,icorr); 
%y = deconvwnr(signal,psf,nsr); 
[y2,x2] = xcorr(signal,psf,127); plot(x2+129,y2,'b.','DisplayName','Cross-Correlation');
y = deconvwnr(signal,psf,nsr); plot(dsp.t,y,'r.','DisplayName','Weiner Deconvolution');


%SUPERSAMPLED SPLINE ------------------------------------------------------
nss = 1000; %number of supersamples
method = 'spline';  
xi = linspace(0,max(dsp.t),nss); 
yi = interp1(dsp.t,y,xi,method);
plot(xi,yi,'r-','color',[1 .7 .7],'Display','Supersample')
title(sprintf('Wiener Deconvolution, %.2g NSR',nsr));

%POINTS ABOVE THRESHOLD ---------------------------------------------------
yi(yi<threshold)=0;
s = diff([0 [0 sign(diff(yi))]]); %second derivative of sign
i = find(s==-2) - 1; %find locations where sign changes from 1 to -1, where s=-2
i = i(i > (nss/dsp.samples) & i < (nss-nss/dsp.samples));
%plot(xi(i),yi(i),'mo') %supersample peaks above threshold

%ERRORS -------------------------------------------------------------------
that = xi(i); %estimated timestamps
dt = that-PE.t(1:ni);
d = minabs(dt,1);

currentitle = get(get(gca,'title'),'string');
title(h(4),sprintf('%s\n%.0f/%.0f PE found, error = %.3f \\pm %.3f ns',currentitle,numel(d),ni,mean(d),std(d)));
hp=fcnplotline(h(2),[nan threshold]); set(hp,'color',[.8 .8 .8],'LineWidth',2,'DisplayName','Threshold');
hp=fcnplotline(h(4),[nan threshold]); set(hp,'color',[.8 .8 .8],'LineWidth',2,'DisplayName','Threshold');
sca(h(4));  title('Weiner Deconvolution');

%PLOT TRUE TIMES ----------------------------------------------------------
fcntight(h,'jointx')
if ~MTCdataflag
    for j=[1 2 4]
        plot(h(j),[1 1]'*PE.t(1:ni)',get(h(j),'ylim')'*ones(1,ni),'-','color',[.5 .5 1])
    end
    xyzlabel(h,'T (ns)','signal (V)')
else
    xyzlabel(h,'T (bins)','signal (bins)')
    delete(h([1 3]));
end
fcnmarkersize(10)






