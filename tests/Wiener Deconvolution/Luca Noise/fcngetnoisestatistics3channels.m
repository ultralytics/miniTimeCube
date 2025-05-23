% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [ output_args ] = fcngetnoisestatistics3channels()
clc; close all
h=fig(2,3,1,1);

% %EXAMPLE DATA -------------------------------------------------------------
% t=linspace(4,4.5,300);
% y=sin(2*pi*10*t)+2*sin(2*pi*50*t)+1.1; %10Hz and 5Hz
% y=y + randn(size(y))*1;

load Y.mat
c = {'r.-','g.-','b.-'};
for i=1:3
    %LOAD DATA ----------------------------------------------------------------
%     load Luca3noise.mat
%     Y(:,1) = x000;
%     Y(:,2) = x212;
%     Y(:,3) = x314;
    y = Y(:,i);
    y = diff(y);
    %y=y(abs(y)<20);
    y = y(y~=0);
    t = (1:length(y));%*372E-12;
    

    %DATA ANALYSIS ------------------------------------------------------------
    dt = t(2)-t(1);
    Fs = 1/dt;
    N = length(y);
    
    %SIGNAL -------------------------------------------------------------------
    sca(h(1))
    plot(t,y,c{i}); xyzlabel('time (bins)','signal (bins)'); title ('signal');
    
    %FIRST DERIVATIVE ---------------------------------------------------------
    sca(h(2))
    plot(t(1:end-1),diff(y),c{i}); xyzlabel('time (bins)','signal (bins)'); title ('derivative');
    
    %AUTOCORRELATION ----------------------------------------------------------
    sca(h(3))
    [acf, lags] = autocorr(y,N-1);
    plot(lags,acf,c{i}); xyzlabel('time lag (bins)','autocorrelation'); title('autocorrelation');
    
    %FIT EXPONENTIAL
    j = find(acf<0,1,'first');  j=60; yi = acf(1:j-1)';
    xi = (0:numel(yi)-1)*.372;
    %fitresult = expFit(xi, yi);
    %-1/fitresult.b
    %plot(xi,yi,'r')
    
    %FFT ----------------------------------------------------------------------
    sca(h(4))
    [y1,f1,w1] = fcnfft(y,t);
    plot(f1,y1,c{i}); xyzlabel('frequency (Hz)','amplitude'); title('FFT');
    
    %PERIODOGRAM --------------------------------------------------------------
    sca(h(5))
    [y1, f1] = periodogram(y-mean(y),[],[],Fs);  % Uses default window, overlap & NFFT.
    plot(f1,y1,c{i}); xyzlabel('frequency (Hz)','power (dB)'); title('Periodogram PSD Estimate');
    
    %WELCH --------------------------------------------------------------------
    sca(h(6))
    [y1, f1] = pwelch(y-mean(y),[],[],[],Fs);  % Uses default window, overlap & NFFT.
    plot(f1,y1,c{i}); xyzlabel('frequency (Hz)','power (dB)'); title('Welch PSD Estimate');
end
sca(h(1))
legend('Channel 1','Channel 2','Channel 3')
fcnfontsize(10)
fcnbox(h,'on')

% %PSD USING FFT ------------------------------------------------------------
% sca(h(5))
% xdft = fft(y-mean(y));
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)).*abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(y):Fs/2;
% plot(freq,10*log10(psdx),'b.-');
% title('FFT PSD Estimate');
% xlabel('Frequency (Hz)'); ylabel('Power / Frequency (dB / Hz)');

% %MULTITAPER --------------------------------------------------------------
% sca(h(7))
% [y1, f1] = pmtm(y-mean(y),[],[],Fs);  % Uses default window, overlap & NFFT. 
% plot(f1,y1,'b.-'); xyzlabel('frequency (Hz)','power (dB)'); title('Periodogram PSD Estimate');
% 
% %YULE WALKER AR --------------------------------------------------------------
% sca(h(8))
% [y1, f1] = peig(y-mean(y),10,N,Fs);  % Uses default window, overlap & NFFT. 
% plot(f1,y1,'b.-'); xyzlabel('frequency (Hz)','power (dB)'); title('Periodogram PSD Estimate');
% 
% %BURG --------------------------------------------------------------
% sca(h(9))
% [y1, f1] = pmusic(y-mean(y),10,N,Fs);  % Uses default window, overlap & NFFT. 
% plot(f1,y1,'b.-'); xyzlabel('frequency (Hz)','power (dB)'); title('Periodogram PSD Estimate');

fcngrid(h,'on')
fcntight(h(4:end),'jointx')
fcntight(h([1 2 3]),'x')
end

function [fitresult, gof] = expFit(xi, yi)
%CREATEFIT(XI,YI)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xi
%      Y Output: yi
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Feb-2014 19:47:42


% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xi, yi );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'Off';
opts.StartPoint = [0.229929895258881 -0.781506121960619];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'yi vs. xi', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'xi' );
ylabel( 'yi' );
grid off
end

