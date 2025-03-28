% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = fiberNeutronFitter(input,output,flags,PE,G1,handles,plotflag)
results.protonprotonflag = [0 0];
results.minangleerror = 0;
results.true = zeros(1,29);
results.xhat = zeros(1,29);
A=output(1);  amplitudes=max(A.D,[],1);  ica = input.cube.all;  Lr=input.cube.Lr;  v128=1:128;  v129=v128+128;

if flags.status.lappd
    captureFraction = 0.25*2; %for air-clad fibers
    PE2E = 1/(input.Material(5).yield*input.cube.QEmean*captureFraction); %PE to MeV
    pid=find(amplitudes(:)>.02 & amplitudes(ica.opposingPixelLAPPD)>.02 & amplitudes(ica.matchedPixelLAPPD)>.02 & amplitudes(ica.opposingPixelLAPPD(ica.matchedPixelLAPPD))>.02)';
    pxyz = ica.xyzLAPPD;
    pid = unique([pid, ica.opposingPixelLAPPD(pid)']);
    pida = pid(ica.fiLAPPD(pid)==5 & ica.leftsideLAPPD(pid));  n=numel(pida);  if n<2; fprintf('<2 scatters found.\n'); return; end %upper pixels
    pidb = ica.matchedPixelLAPPD(pida)';
    pidc = ica.opposingPixelLAPPD(pida)'; %lower pixels
    pidd = ica.matchedPixelLAPPD(pidc)';
    if isempty(input.NN)
        switch input.cube.prettyvolume
            case '146L'
                input.NN = load('NN.TS FiberTrain - 146L LAPPD.80cm 50pc 2pix EJ-254 100k.mat');
            case '7L'
                input.NN = load('NN.TS FiberTrain - 7L LAPPD.20cm 50pc 2pix EJ-254 100k.mat');
        end
    end
    
    %GET VERTICES
    x = zeros(n,5);  t=zeros(n,4);
    [~, ~, ~, t(:,1)] = fcnpulsewidth(A.D(:,pida),.5,[0 1000],[0 1000],1E6,'fraction');
    [~, ~, ~, t(:,2)] = fcnpulsewidth(A.D(:,pidb),.5,[0 1000],[0 1000],1E6,'fraction');
    [~, ~, ~, t(:,3)] = fcnpulsewidth(A.D(:,pidc),.5,[0 1000],[0 1000],1E6,'fraction');
    [~, ~, ~, t(:,4)] = fcnpulsewidth(A.D(:,pidd),.5,[0 1000],[0 1000],1E6,'fraction');
    for i=1:n
        ia=round(min(t(i,:)));  ia=floorandceil(ia-10,0,127);  Vb=fcncol(A.D(:,[pida(i) pidb(i) pidc(i) pidd(i)]));   Va=Vb(ia+[v128 v128+256 v128+256*2 v128+256*3]);  Vb=Vb(1:2:end);
        tbias=[0 0 ia*.2 0];
        
        z=input.NN.neta(Va)+tbias';
        %z=input.NN.netg([input.NN.neta(Va)'+tbias, input.NN.netb(Va)', input.NN.netc(Va)'+tbias(3), input.NN.netd(Vb)', input.NN.nete(Vb)', input.NN.netf(Va)', sum(Vb)]');
        
        x(i,:) = [pxyz(pida(i),1), z(1:3)', z(4)*PE2E]; %[x y z t e]
    end
    x(:,1)=floorandceil(x(:,1),-Lr(1),Lr(1));
    x(:,2)=floorandceil(x(:,2),-Lr(2),Lr(2));
    x(:,3)=floorandceil(x(:,3),-Lr(3),Lr(3));
    x(:,5)=max(x(:,5),0);
    [~, si] = sortrows(x,-5);
    
    %MERGE NEIGHBORS
    anodetol = 15; %norm(input.NN.sigmas(1)*[s s]); %15
    fibertol = 60; %norm(input.NN.sigmas(2)*[s s]); %60
    timetol  = inf; %norm(input.NN.sigmas(3)*[s s]); %inf
    removed  = false(1,n);
    for i=si'
        if removed(i); continue; end
        dx = abs(x(i,:)-x);  j=dx(:,1)<7 & dx(:,2)<anodetol & dx(:,3)<fibertol & dx(:,4)<timetol & rangec(dx(:,1:3))>.01;
        
        if any(j)
            wi=x(i,5);  wj=x(j,5)';  sw=wi+sum(wj);  x(i,1:4) = (wi*x(i,1:4) + wj*x(j,1:4))./sw;  x(i,5) = x(i,5) + sum(x(j,5));  %weighted mean
            removed(j) = true;
        end
    end
    
else
    pid=find(amplitudes>.004 & amplitudes(ica.opposingPixel)>.004)';
    pxyz = ica.xyz;
    pid = unique([pid, ica.opposingPixel(pid)]);
    pida = pid(ica.fi(pid)==5);  n=numel(pida);  if n<2; fprintf('<2 scatters found.\n'); return; end %upper pixels
    pidb = ica.opposingPixel(pida); %lower pixels
    %if isempty(input.NN);  input.NN.fiberfit=load('NN.NORMALIZED.TS FiberTrain 2mm 20cm 3MeV - 2L SENSL J60035 5V 50pc 512pix EJ-254 100k.mat'); end
    
%     %1. GET VERTICES W/ TOF
%     x = zeros(n,5);  D = input.cube.Lr(3)*2;  v=input.Material(5).mu(6);  v=70; %v down fiber (mm/ns)
%     for i=1:n
%         ta = A.t(A.pid==pida(i));  na = A.N(pida(i));
%         tb = A.t(A.pid==pidb(i));  nb = A.N(pidb(i));  if na<3 || nb<3 || isempty(ta) || isempty(tb); x(i,:)=nan; continue; end
%     
%         that = (ta+tb-D/v)/2; %(ns)
%         xyzhat = [pxyz(pida(i),1:2), v*ta/2-v*tb/2]; %(mm)
%         x(i,:) = [xyzhat that (na+nb)*PE2E];
%     end
    
    %1. FIND XYZTE
    t=zeros(n,2);  x=zeros(n,5);  s=x;  NN=input.NN.fiberfit;
    [~, ~, ~, t(:,1)] = fcnpulsewidth(A.D(:,pida),.5,[0 1000],[0 1000],1E6,'fraction');
    [~, ~, ~, t(:,2)] = fcnpulsewidth(A.D(:,pidb),.5,[0 1000],[0 1000],1E6,'fraction');
    ia=round(min(t,[],2));  ia=floorandceil(ia-30,0,127);  Vb=[A.D(1:256,pida)' A.D(1:256,pidb)'];
    
    Va=zeros(n,256);  for i=1:n; Va(i,:)=Vb(i,ia(i)+[v128 v128+256]); end
    Va(:,v128)=Va(:,v128)./max(Va(:,v128),[],2);
    Va(:,v129)=Va(:,v129)./max(Va(:,v129),[],2);
    xa=NN.net{1}(Va')';  xa(:,2)=xa(:,2)+ia*.2;
    xb=NN.net{2}([NNwaveformStats(Vb(:,1:256),2) NNwaveformStats(Vb(:,257:end),2) xa(:,1)]')';
    x = [pxyz(pida,1:2), xa, xb];
    x(:,3)=floorandceil(x(:,3),-Lr(3),Lr(3));  x(:,5)=max(x(:,5),0);
    if ~flags.status.MC
        [~, si] = sortrows(x,4);
        popoutsubplot(handles.GUI.axes1, fig(1,1,2))
        plot3(x(si,1),x(si,2),x(si,3),'k.','markersize',25,'linewidth',2); for i=1:n; j=si(i); text(x(j,1),x(j,2),x(j,3),sprintf(' %g\n',i),'FontSize',14); end
    end
    [~, si] = sortrows(x,-5);
    
    
    %2. FIND UNCERTAINTY
    E=x(:,5);
    s(:,1:2)=sqrt(1/12)*input.Volume.L(1)*2; %uniform distribution std in fiber x,y
    s(:,3:5)=[NN.S{3,1}.s(E), NN.S{3,2}.s(E), NN.S{3,3}.s(E)];  %[ZTE] Uncertanties (sigmas)
    
    if strcmp(G1.pname{1},'mu-'); results=fiberMuonFitter(input,output,flags,PE,G1,handles,x,s,plotflag); return; end
    
    
    %3. MERGE NEIGHBORS
    removed=false(n,1);  assignment=zeros(1,n);
    xr=bsxrangec(x(:,1:2),x(:,1:2));  xr(xr<.01)=inf;  %xr=min(xr); %min range to nearest pixel
    diagonal=xr>8 & xr<9;  neighbors=xr<7; ai=0;
    for i=si'
        if removed(i); continue; end
        sigi = sqrt(s(i,:).^2 + s.^2) * 9; %6sigma merging bounds
        dx = abs(x(i,:)-x);  j=dx(:,1)<6.5 & dx(:,2)<6.5 & dx(:,3)<sigi(:,3) & dx(:,4)<sigi(:,4) & ~removed;
        anysharedneighbors = any(neighbors(:,i) & neighbors);  j(diagonal(i,:) & ~anysharedneighbors) = 0; 
        ai=ai+1; assignment(j)=ai;
        
        if sum(j)>1
            wj=1./s(j,1:4).^2;
            x(i,1:4)=sum(wj.*x(j,1:4))./sum(wj);            x(i,5)=sum(x(j,5));  %inverse-variance weighted mean
            s(i,1:4)=sqrt(1./sum(wj));                      s(i,5)=sqrt(sum(s(j,5).^2));
            removed(j)=true;  removed(i)=false;
        end
    end
end
x=x(~removed,:);
s=s(~removed,:);

%3. VERIFICATION AND PLOTTING
PE.triggeredFibers=numel(pida);
results = doubleScatterVerification(results,input,G1,PE,handles,x,s,G1.pname{1},plotflag);

