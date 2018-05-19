function results = doubleScatterVerification(results,input,G1,PE,handles,x,s,particle,plotflag)
MTCflag = ischecked(handles.GUI.realdataflag);
fiberCaptureFraction=.5;  E2PE=input.Material(5).yield*input.cube.QEmean*fiberCaptureFraction; %MeV to PE
if ~isfield(PE,'triggeredFibers'); PE.triggeredFibers=[]; end

x = sortrows(x,4); %sort by time
nRecoils=size(x,1);  if nRecoils<2; fprintf('<2 scatters found.\n'); return; end
%if abs(x(2,4)-x(1,4))<sqrt(s(2,4).^2+s(1,4).^2)*0.6 && ~MTCflag; fprintf('dt within 75pc 1sigma uncertainty'); return; end

%PARTICLE ID
Nhat = x(1:2,5)'*E2PE; %Sum Recoil 1 & 2 PEs
vhat = rangec(x(2:end,1:3),x(1:end-1,1:3))./(x(2:end,4)-x(1:end-1,4)); %mm/ns
E0unc = sum(x(:,5)); % uncorrected sum! 
aStats = [E0unc, PE.triggeredFibers, mean(diff(x(:,4))), var(x,[],1), mean(vhat), var(vhat)];
X = [x(1,1:4) x(2,1:4) Nhat nRecoils vhat(1) aStats]; %X=[X(i,[1:4  6:9  11:12 14:15])  X(i,18:end)];
P=[0 0];  Ptrue = double([strcmp(particle,'neutron') strcmp(particle,'gamma')]);
%if input.tsv==2
    p=.95; %p=probability threshold
    %P = input.NN.fiberpid.dnet(X')'; %P=particleProbability
    P = Ptrue;
    if P(1)>p
        particle='neutron';    fprintf('%.4f neutron PID probability.\n',P(1));
    elseif P(2)>p
        particle='gamma';      fprintf('%.4f gamma PID probability.\n',P(2));
    else
        fprintf('Ambiguous PID, [%.3f %.3f] probability.\n',P(1),P(2)); return;
    end
%end

switch particle
    case 'neutron'
        dEhat = E2neutronE(x(:,5)); %proton recoil MTC calibration curve
        E1hat = fcnneutronke(x); %neutron ke after first bounce (MeV) (different from gamma value)
        E0hat = E1hat + dEhat(1); %original neutron energy, E0hat = sum(dEhat)./(1-.5^nRecoils);
        anglehat = fcnneutrontheta(E0hat,dEhat(1)); %half angle
        color='r';
    case 'gamma'
        %COMPTON SCATTER EQUATION  https://en.wikipedia.org/wiki/Compton_scattering -------
        % c = 299792458; %m/s, speed of light
        % h = 4.13566751691E-15; %eV*s %planck constant
        me = 0.510998910; %electron rest mass (MeV) or (MeV/c^2)
        % syms h c E2 E1 me theta;  solve(h*c/E2 - h*c/E1 - h/me/c*(1-cos(theta)),theta)
        dEhat=x(:,5);          
        %E0unc = sum(dEhat)./(1-.5^nRecoils);  %%assume it loses half its energy every recoil

        if nRecoils>200 && input.tsv==2
            %syms E0 dE me theta;  solve(1/(E0-dE) - 1/E0 - 1/me*(1-cos(theta)),E0)
            dE = x(2,5); %second recoil
            theta = fcnangle(x(2,1:3)-x(1,1:3),x(3,1:3)-x(2,1:3)); %(rad)
            %E1 = ((-dE*(cos(theta) - 1)*(dE + 4*me - dE*cos(theta)))^(1/2) - dE + dE*cos(theta))/(2*(cos(theta) - 1));
            E1 = -(dE + (-dE*(cos(theta) - 1)*(dE + 4*me - dE*cos(theta)))^(1/2) - dE*cos(theta))/(2*(cos(theta) - 1));
            E0hat = gammaE0(input,x,[E1+x(1,5) s(2,5)*2]);
        else
            E0hat = gammaE0(input,x);
        end

        E1=E0hat;  E2=E0hat-dEhat(1);  anglehat=acosd( me*(1./E1-1./E2)+1 );
        ehat=acosd( me*(1./E0hat-1./(E0hat-dEhat(1)))+1 );
        color='b';
    otherwise %pointsource

end
if imag(anglehat)~=0; fprintf('imaginary angle cone.\n'); anglehat=real(anglehat); end

if ~MTCflag
    nv = find(G1.tid==1);  if numel(nv)<3; fprintf('insufficient GEANT TIDs.\n'); return; end %neutron rows
    nRecoilsTrue = sum(G1.p2inside(nv)); %total recoils inside detector
    p0 = G1.p1(nv(1),:); %(mm) where it enters the detector
    p1 = G1.p1(nv(2),:); %first bounce
    p2 = G1.p2(nv(2),:); %second bounce
    
    t1 = G1.t2(nv(1)); %(ns) first bounce
    t2 = G1.t2(nv(2)); %second bounce
    v = rangec(p1,p2)./(t2-t1); %mm/ns
    
    En = G1.ke1(nv);
    e1 = En(1)-En(2); %(MeV) first bounce dE
    e2 = En(2)-En(3); %(MeV) second bounce dE
    
    vec1 = p1 - p0; %neutron vec before bounce
    vec2 = p2 - p1; %neutron vec after bounce
    angle = fcnangle(vec2,vec1)*r2d; %rad
    %angleFromEn=acosd( me*(1./En(1)-1./En(2))+1 );
    N = [sum(PE.Gptid==2) sum(PE.Gptid==3)]; %photons observed from first and second bounce
    
    results.protonprotonflag = [any(G1.pid(G1.tid==2)==2212) any(G1.pid(G1.tid==3)==2212)];
    results.minangleerror = fcnangle(vec1,x(2,1:3)-x(1,1:3))*r2d - anglehat;
    
    Ptrue=0;
    results.true = [reshape([p1 t1 e1; p2 t2 e2]', [1 10]),       N,    angle, nRecoilsTrue,      v,                     0, En(1), En(1),                  0,                  0,  zeros(1,5),          0,        0   Ptrue  ];
else
    En=0; angle=0;
    results.protonprotonflag = [false false];
    results.minangleerror = 0;
    results.true = zeros(1,28);
end
results.xhat = [reshape([x(1:2,1:4) dEhat(1:2)]', [1 10]), Nhat, anglehat, nRecoils,    vhat(1), results.minangleerror, E0hat, E0unc, PE.triggeredFibers, mean(diff(x(:,4))), var(x,[],1), mean(vhat), var(vhat)  P      ];

if plotflag
    fprintf('[%.3f fit E0, %.3f true],  [%.3f fit angle, %.3f true]\n',E0hat,En(1),anglehat,angle)
    closeallexcept(handles.GUI.figure1);
    popoutsubplot(handles.GUI.axes1, fig(1,1,2))
    p1hat = x(1,1:3);
    p2hat = x(2,1:3);
    [cx,cy,cz] = fcnplotcone(p1hat,p2hat,anglehat);
    plot3(x(:,1),x(:,2),x(:,3),[color '.-'],'markersize',25,'linewidth',2)
    plot3(x(:,1),x(:,2),x(:,3),'g.-','markersize',25,'linewidth',2); for i=1:nRecoils; text(x(i,1),x(i,2),x(i,3),sprintf(' %g\n',i),'FontSize',14,'Color','g'); end

    surf(cx,cy,cz,'facecolor',color,'edgecolor','none','facealpha',.3); axis equal vis3d; %camlight headlight;
    
    fprintf('\n     1x     2y     3z     4t     5E     6x     7y     8z     9t    10E    11N    12N  13angle  14np   15v    16aerror 17E0  18E0unc 19Trig  20-24Var(x)\n');
    fprintf(' %6.1f',results.xhat); fprintf('\n'); fprintf(' %6.1f',results.true); fprintf('\n');
end
end

function E0hat=gammaE0(input,x0,a)
%https://en.wikipedia.org/wiki/Klein?Nishina_formula
%a = apriori [mu sigma]
np=numel(x0)/5;  zf=cell(np+1,1);  xf=linspace(0,10,10000)';
for i=1:np
    y=double( input.optimizer.gamma.F({xf,x0(i,5)}) );
    if i>1
        y=interp1c(xf+x0(1,5),y,xf);
        %y=interp1(xf+x0(1,5),y,xf);
    end
    zf{i}=y;
end
if nargin>2;  zf{np+1}=pdf('norm',xf,a(1),a(2));  end
zf=cat(2,zf{:}); p=prod(zf,2);

%fig(2,1,'19cm'); sca; plot(xf,zf); xyzlabel('E_\gamma (MeV)'); sca; plot(xf,p);  xyzlabel('E_\gamma (MeV)'); fcntight('equal'); fcnlinewidth(2)
[~,i]=max(p);  E0hat=xf(i);
end


function x=fcnNLS(x0,s)
sx=size(x0');  x0=fcncol(x0')';  s=fcncol(s')';   z=x0;

options = optimoptions('fminunc','Display','off','MaxIter',400,'Algorithm','quasi-newton'); % run interior-point algorithm

x = fminunc(@(x) fcnL(z,s,x), x0, options);

x=reshape(x(:),sx)';

% for i=1:3
%     J = eye(n);
%     zhat = 
%     x = x + (J'*R*J)^-1*J'*R*(z-zhat)*damping;
% end
end



function fx = fcnL(mu,sigma,x)
%X = [xyzte, xyzte ...];
n=numel(x);  np=n/5;
iv      = (0:np-1)'*5;
posi    = (1:3)+iv;  ti=4+iv;
bi=2:np;  ai=1:np-1;
dx      = rangec(x(posi(bi,:)),x(posi(ai,:)));
dt      = x(ti(bi)) - x(ti(ai));
vel     = dx(:)./dt(:);

%y = log( exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma) );
y = -0.5 * ((x - mu)./sigma).^2 - log(sqrt(2*pi) .* sigma);

sigmac = 100;
yc = -0.5 * ((vel - 299.792458)./sigmac).^2 - log(sqrt(2*pi) .* sigmac);

fx = -sum(y)-sum(yc);
end
