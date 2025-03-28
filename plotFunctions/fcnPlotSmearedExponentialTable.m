% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function fcnPlotSmearedExponentialTable(input)
fig(1,1,'15cm');
a = input.table.smearExp;
plot(a.x, a.pdfsmeared, '.-');
xyzlabel('T (ns)','pdf','',sprintf('%s scintillator response\nt_r = 0.85 ns, t_f = 1.51 ns, \\sigma = 0.10 ns',input.Material(1).name))
fcntight
h=gca; h.XLim(2)=15;

%[~,amplitude,w,t,int]=fcnpulsewidth(a.pdfsmeared,.5,[0 1000],[0 1000],1E6,'fraction')
