function [] = fcnKleinNishina()
clc; clear all;  close all
me = 0.510998910; %electron rest mass (MeV) or (MeV/c^2)
re=1; a=1;
%theta=linspace(0,pi,1000);

n=2000; ev=linspace(0,10,n); X=zeros(n);
for i=1:n
    E1 = ev(i);  % initial gamma enerrgy
    E2 = linspace(0,E1,n);  % gamma energy after scatter
    theta=real( acos( me*(1./E1-1./E2)+1 ) );  % gamma scatter angle
    
    j=theta==pi;
    
    %syms E2 E1 me theta;  solve(1/E2 - 1/E1 - 1/me*(1-cos(theta)),E1)
    %ratio = 1./(1 + (E1./me)*(1-cos(theta))); %E2./E1
    
    %dsdTh = re^2/2 * (E2./E1).^2.*(E2./E1+E1./E2 - sin(theta).^2); %https://en.wikipedia.org/wiki/Klein?Nishina_formula
    %y2 = a^2*re^2*ratio.^2.*(ratio+1./ratio-1+cos(theta).^2)/2;
    %dsdT = dsdTh.*(2*pi./(E1-(E1-E2)).^2);
    dsdT = (re^2*pi*(E1./E2 + E2./E1 - sin(real(acos(me*(1./E1 - 1./E2) + 1))).^2))./E1.^2; %simplified version

    dsdT(j)=0; theta(j)=0;
    %sca; plot(E1-E2,theta*r2d,'.'); sca; plot((E1-E2),dsdT,'.'); %plot(E2,y2,'.')
    
    if E1>0
        X(i,:)=interp1((E1-E2),dsdT,ev);
    end
end
X=single(X);

X=X./nansum(X,2);
F=griddedInterpolant({ev,ev},X,'linear','none');
fig(2,1,'19cm') 
for a=[.511 1 1.5 2]
    plot(ev,F(zeros(size(ev))+a,ev),'DisplayName',sprintf('%.3f MeV \\gamma',a));
end
xyzlabel('E_{e^-} (MeV)'); legend show; fcnlinewidth(2); fcnfontsize(12); title('Estimating e^- Energy From \gamma')

sca;
X=X./nanmean(X,2);
F=griddedInterpolant({ev,ev},X,'linear','none');
for a=[1.81 .07]
    plot(ev,F(ev,zeros(size(ev))+a),'DisplayName',sprintf('%.3f MeV e-',a));
end
xyzlabel('E_\gamma (MeV)'); legend show; fcnlinewidth(2); fcnfontsize(12); title('Estimating \gamma Energy From e^-')

fig; imagesc(ev,ev,X'); fcntight('c'); xyzlabel('E_\gamma (MeV)','E_{e^-} (MeV)'); fcntight('xy csigma'); xyzlabel('E_\gamma (MeV)','E_{e^-} (MeV)')
%save gammaKleinNishinaF.mat F -v6
end



function E0hat=plotKN()
%https://en.wikipedia.org/wiki/Klein?Nishina_formula
E=[1173.2 1332.5]/1E3;
load gammaKleinNishinaF.mat
np=numel(E);  zf=cell(np+1,1);  xf=linspace(0,max(E)*1.1,10000)';
for i=1:np
    zf{i}=double( F({E(i),xf}) );
end
zf=cat(1,zf{:}); p=sum(zf,1);

fig(2,1,'19cm'); sca; plot(xf,zf); xyzlabel('dE_\gamma (MeV)'); sca; plot(xf,p);  xyzlabel('dE_\gamma (MeV)'); fcntight('equal'); fcnlinewidth(2)
end
