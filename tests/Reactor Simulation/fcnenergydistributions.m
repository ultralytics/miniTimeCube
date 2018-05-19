function spectra = fcnenergydistributions()
plotflag=0;
ls = cell(4,1);
x = linspace(.001,11,1000); dx = x(2)-x(1);

%ANTINEUTRINOS ------------------------------------------------------------
s = fcnspec1s(x, .005, [], 1, 1, 1);
y = s.s/sum(s.s(~isnan(s.s)))/dx;
if plotflag; figure; loglog(x,y,'k'); hold on; end
ls{1} = sprintf('Antineutrinos: 5m range pdf (Fogli + Daya Bay 2012)');
spectra(1).pdfx = x;  spectra(1).pdfy = y;

%REACTOR NEUTRONS ---------------------------------------------------------
k = 1.5; 
theta = 0.64;
y = pdf('gamma',x,k,theta);
if plotflag; loglog(x,y,'r'); end
p = sum(y(1:199))/sum(y)*100; %percent under 2MeV
ls{2} = sprintf('Reactor Neutrons: gamma pdf, k=%3g, \\theta=%3g (%.1f%%<2MeV)',k,theta,p);
spectra(2).pdfx = x;  spectra(2).pdfy = y;

%ATMOSPHERIC NEUTRONS -----------------------------------------------------
load Goldhagenneutron.mat ; %y = (cm^-2 s^-1)
y=y./x; %y = (cm^-2 s^-1 MeV^-1)
if plotflag; loglog(x,y,'g'); end
ls{3} = sprintf('Atmospheric Neutrons: (Goldhagen 2004)');
spectra(3).pdfx = x';  spectra(3).pdfy = y';

% pdfx = x;  pdfy = yreactor*MWth;
% cdfy = cumsum(pdfy.*[0; diff(pdfx)]);
% y = fcnrandcdf(cdfy,pdfx,1E4,'linear');

%REACTOR GAMMAS -----------------------------------------------------------
% lambda = 6.5; 
% k = 4;
% y = pdf('wbl',x,lambda,k);
% if plotflag; semilogx(x,y,'g'); end
% ls{4} = sprintf('Gammas: weibull pdf, \\lambda=%3g, k=%3g',lambda,k);
% spectra(4).pdfx = x;  spectra(4).pdfy = y;
MWth=20;
load gamma.mat
y = (ybackground*0+yreactor*MWth)';
if plotflag; loglog(x,y,'b'); end
ls{4} = sprintf('Reactor Gammas: Pieter Mumm measurements, NIST, %.0fMWth Reactor',MWth);
spectra(4).pdfx = x;  spectra(4).pdfy = y;

%ATMOSPHERIC GAMMAS -------------------------------------------------------
%y = (ybackground+yreactor*MWth*0)';
load('hist_gamma_mTC_9-25-15.mat')
x = linspace(.001,101,8192)'; dx = x(2)-x(1);
y = interp1(bin,hist_gamma_flux./bin,x)';

if plotflag; loglog(x,y,'c'); end
ls{5} = sprintf('Atmospheric Gammas: Pieter Mumm measurements, NIST');
spectra(5).pdfx = x;  spectra(5).pdfy = y;

% pdfx = x;  pdfy = yreactor*MWth + ybackground;
% cdfy = cumsum(pdfy.*[0; diff(pdfx)]);
% y = fcnrandcdf(cdfy,pdfx,1E4,'cubic');


%MUONS --------------------------------------------------------------------
%Source http://muradio.pd.infn.it/bibliografia_files/rastin_JoPG10_1984_1609.pdf
ls{6} = sprintf('Muons: Sea Level pdf (Rastin 1984)');

E_GeV=[0.35
    0.40
    0.50
    0.60
    0.70
    0.80
    0.90
    1.00
    1.50
    2.00
    3.00
    4.00
    5.00
    6.00
    7.00
    8.00
    9.00
    10.00
    15.00
    20.00
    25.00
    30.00
    40.00
    50.00
    60.00
    70.00
    80.00
    90.00
    100.00
    150.00
    200.00
    250.00
    300.00
    400.00
    500.00
    600.00
    700.00
    800.00
    900.00
    1000.00
    1500.00
    2000.00
    2500.00
    3000.00]';

%(cm-1 sr-1 s-1 (GeV/c)-1)
Muon_differential_intensity=[2.85E-3
                            2.90E-3
                            2.94E-3
                            2.92E-3
                            2.87E-3
                            2.80E-3
                            2.71E-3
                            2.62E-3
                            2.12E-3
                            1.69E-3
                            1.10E-3
                            7.40E-4
                            5.17E-4
                            3.75E-4
                            2.80E-4
                            2.16E-4
                            1.69E-4
                            1.35E-4
                            5.28E-5
                            2.58E-5
                            1.45E-5
                            8.69E-6
                            3.90E-6
                            2.11E-6
                            1.26E-6
                            8.03E-7
                            5.42E-7
                            3.81E-7
                            2.77E-7
                            7.85E-8
                            3.12E-8
                            1.50E-8
                            8.20E-9
                            3.11E-9
                            1.45E-9
                            7.75E-10
                            4.55E-10
                            2.86E-10
                            1.89E-10
                            1.31E-10
                            3.14E-11
                            1.13E-11
                            5.11E-12
                            2.67E-12]';

% pdfx = log10(E_GeV);  pdfy = Muon_differential_intensity;  fig; plot(cdfy,pdfx,'.-')
% cdfy = cumsum(pdfy.*[0 diff(pdfx)]);
% y = 1000 * 10.^fcnrandcdf(cdfy,pdfx,1E4,'pchip');  %(MeV)
% fig; hist(log10(y),30);

%Normalize into PDF
Muon_intensity=Muon_differential_intensity;
average_muon_energy_at_sea_level=sum(E_GeV.*Muon_intensity)/sum(Muon_intensity);
fig; plot(log(E_GeV),Muon_intensity,'b-o')
spectra(6).pdfx = E_GeV*1E3;  spectra(6).pdfy = Muon_intensity;


for i=1:6
    spectra(i).name = ls{i};
end

if plotflag
    loglog(E_GeV*1000,Muon_intensity,'m');
    
    fcnlinewidth(2)
    fcnfontsize(8)
    legend(ls{1:6},'Location','SouthWest')
    xlabel('particle energy (MeV)')
    ylabel('pdf')
    title('Particle Energy Spectra')
    grid on
end
end