% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = fcnantineutrino(input,output,handles,G1,plotflag)
zv=zeros(1,10); results.xhat=zv; results.true=zv;  if plotflag; try closeallexcept(handles.GUI.figure1); catch; end; tic; end

%TRUTH --------------------------------------------------------------------
if ~isempty(G1)
    inu=find(G1.pid==-12,1,'first');
    i=find(G1.pid==-11,1,'first');
    j=find(any(G1.pid==[1000030070, 1000020040, 1000010020, 1000020032, 1000030080],2)); %Li7 or alpha or deuteron or He3 or Li8 
    results.true = [G1.p1(i,:) G1.t1(i) G1.ke1(inu), mean(G1.p1(j,:),1) mean(G1.t1(j)) sum(G1.de(j))];
end
if output(1).Nsum<5 || output(1).Nsum>1E6;  return;  end


%FERMAT ESTIMATORS --------------------------------------------------------
ax=cell(2,1);  ax(:)={nan(1,5)};  bx=ax;
for i=1:numel(output)
    if output(i).Nsum>4  %POISSON ESTIMATOR
        ax{i} = fcnMLpoint(input,output(i).N);
    end
    
    t = output(i).t;  t=t(t~=0);
    if numel(t)>4
        t0=min(fcnsigmarejection(t,2,3));
        k = fcnoptimizerkr(input,output(i));  k.reflections=0;  k.timeflag=1;
        nb = 1;  nx = nb*5;
        x0 = zeros(nb,5);  x0(:,5)=1/nb;  x0(:,4)=t0;  x0(:,1:3)=ax{i}(1:3);  x0=reshape(x0',[1 nx]);
        es = 3; %mm airgap between xhat and detector wall
        Lr = input.cube.Lr - es;
        
        A=[];  B=[];  Aeq=zeros(nb,5);  Aeq(:,5)=1;  Aeq=reshape(Aeq',[1 nx]);  Beq=1;  LB=[];  UB=[];
        if input.cube.shapeID~=3
            LB=repmat([-Lr t0-sqrt(sum(Lr.^2))*k.ci-10  .001],[1 nb]);  UB=repmat([ Lr t0+30 1],[1 nb]);
        end
        %if isnan(fcnfermatpointvectorized(k,output(i).t,x0)); continue; end
        %[bx{i}, fx] = fmincon(@(x) fcnfermatpointvectorized(k,output(i).t(output(i).t~=0),x), x0, A, B, Aeq, Beq, LB, UB, [],input.optimizer.options2);
        bx{i} = patternsearch(@(x) fcnfermatpointvectorized(k,output(i).t(output(i).t~=0),x), x0, A, B ,Aeq, Beq, LB, UB, [],input.optimizer.psoptions3);
        %bx{i} = [0 0 0 t0];
    end
end
results.xhat = [ax{1}(1:3) bx{1}(4) ax{1}(4), ax{2}(1:3) bx{2}(4) ax{2}(4)];

if plotflag==1
    zN=output(1).N;  X=results.xhat(1:5);
    k = fcnoptimizerk(input,zN);
    
    F0      = fcnsolidanglevectorized(k,X(1:3));
    F0      = sum(F0,3); %sum any reflections
    e       = fcnrow(X(5))*k.yield;
    wf      = F0.*e;
    lambda  = zeros(numel(k.fi),1); lambda(k.pid) = wf;

    ha=fig(1,2,2);
    sca;  fcnplotdetectorprojection(input,zN>0,zN);             title('Measured  z')
    sca;  fcnplotdetectorprojection(input,lambda>0,lambda);     title('Reconstructed  z=H(x)')
    [el,az]=fcnelaz(X(1),X(2),X(3));   [rm,cm]=fcnmapprojection(el*r2d,az*r2d,'winkeltripel');  for i=1:2; plot(ha(i),cm,rm,'b.-','markersize',25); plot(ha(i),cm,rm,'b+','markersize',50,'linewidth',1.5); end
    linkprop(ha(1:2),{'CameraPosition','CameraUpVector'}); 
elseif plotflag==2
    ha=fig(2,2,1.5); s={sprintf('  %.2f MeV PROMPT\n  %.4g ns',ax{1}(4),bx{1}(4)),sprintf('  %.1f keV DELAYED\n  %.4g ns',ax{2}(4)*1E3,bx{2}(4))};  rm=zeros(2,1); cm=rm;
    xa=repmat((1:256)',[1 1536]); p{1}=[nan nan nan]; p{2}=p{1};
    for i=1:2
        try
        A=output(i); A.N(A.N==0)=nan;   x=results.xhat((1:5)+(i-1)*5);  A.D(end,:)=nan; 
        sca(i);  plot(xa(:),A.D(:),'color',fcndefaultcolors(i),'linewidth',.6); xyzlabel('T (sample)','V (bins)','',s{i})
        sca(i+2);  fcnplotdetectorprojection(input,A.N>0,A.N);  h=gca;  h.CameraViewAngle=8; h.Position(1)=h.Position(1)-.02;
        p{i}=x;
        end
    end
    fcntight(ha([1 2]),'xyjoint');  ha(2).YAxis.Visible='off';
        
    ax=interp1([0 1],[p{1}(1) p{2}(1)],linspace(0,1,100));
    ay=interp1([0 1],[p{1}(2) p{2}(2)],linspace(0,1,100));
    az=interp1([0 1],[p{1}(3) p{2}(3)],linspace(0,1,100));  [el,az]=fcnelaz(ax,ay,az);   [rl,cl]=fcnmapprojection(el*r2d,az*r2d,'winkeltripel');
    for i=1:2
        try
            sca(ha(i+2));
            plot(cl(1),rl(1),'.','markersize',50,'Color',fcndefaultcolors(1));
            plot(cl(end),rl(end),'.','markersize',50,'Color',fcndefaultcolors(2));  j=[1 numel(rl)];
            text(cl(j(i)),rl(j(i)),s{i},'fontsize',16);
            plot(ha(i+2),cl,rl,'.-','markersize',5,'Color',[.5 .5 .5],'linewidth',1.5); end
        end
    end
end



%ehat = fcnMLenergy(input,x(1:3),output(1).N)


% %POSITRON FIT -------------------------------------------------------------
% r = max(3.5*ehat,3);  %e+ dr/dE = 3.5mm/MeV mean, 3mm min
% xyz = fcnuniformsphere(20,40);
% n=size(xyz,1);  p0 = repmat(x(1:3),[n 1]);
% Xa = [-xyz*r+p0 xyz*r+p0 repmat(x(5),[n 1])];  k = fcnoptimizerkr(input,output(1)); 
% 
% %ITERATE
% Lr=input.cube.Lr-1; A=[];  B=[];  LB=[];  UB=[];  Aeq=[];  Beq=[];  tlb = x(5)-sqrt(sum(Lr.^2))*k.ci;  tub = x(5)+30;
% if input.cube.shapeID~=3; LB=[-Lr -Lr tlb]; UB=[Lr Lr tub]; end %if not cylinder
% %tic; Xb = fmincon(@(x) fermatLine(k,output(1).t,x), Xb, A, B ,Aeq, Beq, LB, UB, [], input.optimizer.psoptions2); toc
% 
% 
% %REGULAR
% k.reflections=0;  k.timeflag=0;  k.cherenkov=0;
% fx = fermatLine(k,output(1).t,Xa);  [~,i]=min(fx);  Xb = Xa(i,:);
% results.xhat(11:13) = Xb(1:3);  [~, results.xhat(14)] = fcnangle(puv,Xb(4:6)-Xb(1:3)); %ct
% 
% 
% %WITH CHERENKOV
% k.reflections=0;  k.timeflag=0;  k.cherenkov=1;
% fx = fermatLine(k,output(1).t,Xa);  [~,i]=min(fx);  Xb = Xa(i,:);
% results.xhat(19:21) = Xb(1:3);  [~, results.xhat(22)] = fcnangle(puv,Xb(4:6)-Xb(1:3)); %ct
% 
% 
% %WITH TIME+CHERENKOV
% k.reflections=0;  k.timeflag=1;  k.cherenkov=1;
% fx = fermatLine(k,output(1).t,Xa);  [~,i]=min(fx);  Xb = Xa(i,:);
% results.xhat(23:25) = Xb(1:3);  [~, results.xhat(26)] = fcnangle(puv,Xb(4:6)-Xb(1:3)); %ct
% %tic; Xb = fmincon(@(x) fermatLine(k,output(1).t,x), Xb, A, B ,Aeq, Beq, LB, UB, [], input.optimizer.psoptions2); toc
% 







function fx = fermatLine(k,zt,xr)
x = reduced2full(xr);
fx = fcnfermatpointvectorized(k,zt,x);
end

function X = reduced2full(Xr)
p1 = Xr(:,1:3);  p2 = Xr(:,4:6);  T = Xr(:,7);  n = 20;  %[xyz xyz t] to [xyztw, xyztw...]
%v=elaz2CC(Xr(:,6:7)); p=Xr(:,1:3);  r=Xr(:,7);  p1=p-v*r/2; p2=p+v*r/2; T=Xr(:,4);  %[xyzt r el az] to [xyztw, xyztw...]
iv = (0:n-1)*5;

v = linspace(0,1,n); %delta points vector
[r,dx] = rangec(p2,p1); %r (mm)
X = zeros(size(Xr,1),5*n);
for i=1:3
    X(:,iv+i) = p1(:,i) - dx(:,i)*v;
end
X(:,iv+4) = T + (r/299.792458)*v; %t (ns), %speed of light 299.79(mm/ns);
X(:,iv+5) = 1/n; %weighti
end