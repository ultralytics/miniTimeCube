% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = gatherTrainingData(input,output,PE,G1)
% mode='trainLAPPDNTC2';
% switch mode
%     case 'trainLAPPDNTC'
%         trainLAPPDNTC
%     case 'TRAIN_FIBER1'
%         trainFiber1
% end

results = [];  A = output(1);  LAPPDflag=regexpi(input.cube.prettyname,'LAPPD');
opposingPixels=input.cube.all.opposingPixel;  fi=input.cube.all.fi;
if PE.n==0; return; end


if LAPPDflag
    B=A.D; [nx, ns, ~, nlappd] = size(B); %number of times, number of strips (30), number sides (2), number of lappd's
    B=reshape(B,[nx ns*2*nlappd]);
    B=fcnSplinePulseReconstruction(B,input.cube.dsp.dynamicrange(2));  %B=reshape(B,[nx ns 2 nlappd]);
    A.D = B;  opposingPixels=input.cube.all.opposingPixelLAPPD;  fi=input.cube.all.fiLAPPD;
end
i=1:256;

%GET CORRECT PIXEL PAIR
pid=mode(PE.pixel(PE.Gpid==-11));
%amplitudes=max(A.D,[],1);  pid=find(amplitudes>.004 & amplitudes(opposingPixels)>.004)';
if numel(pid)==0 || isnan(pid); return; end
pid = unique([pid, opposingPixels(pid)]);
pida = pid(fi(pid)==5); %upper pixels
pidb = opposingPixels(pida); %lower pixels
results.VsumThreshold = [sum(A.D(i,pida),2)',  sum(A.D(i,pidb),2)'];
%n=numel(pida);  if n>=1;  j=randi(n,1);  pida=pida(j);  pidb=pidb(j);  end

if LAPPDflag
    M=input.cube.all.matchedPixelLAPPD;
    if ~input.cube.all.leftsideLAPPD(pida)
        pida=M(pida); pidb=M(pidb);  %pida must be left side!!
    end 
    results.xhat = [A.D(i,pida)' A.D(i,M(pida))' A.D(i,pidb)' A.D(i,M(pidb))'];
    results.true = [G1.p1(1,[2 3]), G1.t1(1), sum(PE.Gptid==1), input.cube.all.xyzLAPPD(pida,1), mean(PE.x(:,1)), min(PE.t(PE.pixel==1)), min(PE.t(PE.pixel==2))]; %[z t E PE] %ONLY FOR POINT SOURCES!!!
else
    results.xhat = [A.D(i,pida)' A.D(i,pidb)'];
    %results.VsumAll = [sum(A.D(i,input.cube.all.fi==5),2)',sum(A.D(i,input.cube.all.fi==6),2)'];  results.xhat = [results.VsumAll, results.VsumThreshold]; %WARNING: SUMS ALL CHANNELS!
    results.true = [G1.p1(1,3), G1.t1(1), A.Ntrue(pida)+A.Ntrue(pidb), G1.de(1)]; %[z t PE E] %ONLY FOR POINT SOURCES!!!
end

end

function trainFiber64_Neutron_Gammadiscrimination()
clc;  clear; close all;  files=uigetfile('*.mat','MultiSelect', 'on');  if ischar(files); files={files}; end;  nf=numel(files); X=[]; T=[]; vertex=[]; %#ok<*NASGU>
for fi=1:nf
    fprintf('Loading ''%s''... ',files{fi});  load(files{fi});  fprintf('Done.\n');
    X=[X; MC.xhat(:,:,1)];  T=[T; MC.xtrue(:,:,1)];  vertex=[vertex;  MC.vertex(:,:,1)]; 
end
clear tsv input
n=size(MC.xhat,1);  T=zeros(size(X,1),2); T(1:n,1)=1; T(n+1:end,2)=1;  n=size(X,1); 

i=~all(X==0,2) & ~all(T==0,2) & ~any(isnan(X),2) & ~any(isnan(T),2);  n=sum(i);
X=[X(i,[1:4  6:9  11:12 14:15])  X(i,18:end-2)];
T=T(i,:);  vertex=vertex(i,:);

[neta, perf, sigmas]=NNtrain(X,T,[10],3); plotconfusion(T',neta(X')); %#ok<*ASGLU>
%P=polyfitrobust(I,T,20);

softnet = trainSoftmaxLayer(neta(X'),T','LossFunction','crossentropy');
dnet = stack(neta,softnet);  dnet.performFcn='crossentropy';
dnet.trainParam.epochs=19000;  dnet.trainParam.time=60*60*4;  dnet.trainParam.max_fail=250;
dnet.divideFcn='dividerand';  dnet.divideParam.trainRatio=70/100;  dnet.divideParam.valRatio=15/100;  dnet.divideParam.testRatio=15/100;
dnet = train(dnet,X',T');  plotconfusion(T',dnet(X')); %#ok<*ASGLU>

e=dnet(X')'-T; i=max(dnet(X')',[],2)>.99; [mean(i(T(:,1)==1))  mean3(i(T(:,2)==1))]
hf=figure; hcp=plotconfusion(T(i,:)',dnet(X(i,:)')); %#ok<*ASGLU>

ha=fig(2,3,'19x28.5cm');
[name, c]=fcnpid2name(2112); cn=c{1};  [name, c] = fcnpid2name(22); cg=c{1};
sca; j=T(:,1)==1; jprob=dnet(X(i & j,:)');  plot(vertex(i & j,5),jprob(1,:)','.','Color',cn); xyzlabel('E_n (MeV)','n Probability')
sca; j=T(:,2)==1; jprob=dnet(X(i & j,:)');  plot(vertex(i & j,5),jprob(2,:)','.','Color',cg); xyzlabel('E_\gamma (MeV)','\gamma Probability')
sca; fcncopyaxes(hf.Children(end),gca); close(hf);
sca; j=T(:,1)==1; [y,x]=histcounts(vertex(j,5),50); y2=histcounts(vertex(j & i,5),x); bar(x(1:end-1)/2+x(2:end)/2,y2./y,1,'FaceColor',cn);  xyzlabel('E_n (MeV)','Acceptance')
sca; j=T(:,2)==1; [y,x]=histcounts(vertex(j,5),50); y2=histcounts(vertex(j & i,5),x); bar(x(1:end-1)/2+x(2:end)/2,y2./y,1,'FaceColor',cg);  xyzlabel('E_\gamma (MeV)','Acceptance')
fcnfontsize(14); delete(sca); fcnmarkersize(6)

%s=sprintf('NN.PATTERN2.%s',files{1});  fprintf('Saving ''%s''... ',s);  save('-v6',s,'dnet');  fprintf('Done.\n\n');
end



function trainLAPPDNTC()
clc;  clear; close all;  files=uigetfile('*.mat','MultiSelect', 'on');  if ischar(files); files={files}; end;  nf=numel(files);
for fi=1:nf
    fprintf('Loading ''%s''... ',files{fi});  load(files{fi});  fprintf('Done.\n');
    inputs=MC.xhat;  targets=MC.xtrue;  clear tsv input
    
    vi{1}=1:256; vi{2}=vi{1}+256; vi{3}=vi{2}+256; vi{4}=vi{3}+256;  amplitudeThreshold=max(inputs(:,vi{1}),[],2)>.02 & max(inputs(:,vi{2}),[],2)>.02 & max(inputs(:,vi{3}),[],2)>.02 & max(inputs(:,vi{4}),[],2)>.02;
    i=find(~all(inputs==0,2) & ~all(targets==0,2) & ~any(isnan(inputs),2) & amplitudeThreshold & MC.xtrue(:,4)>4);  i=i(1:end); inputs=inputs(i,:);  targets=targets(i,:);  n=size(inputs,1);
    
    pmt.stripdx = 6.9088;  pmt.stripx = midspace(-13,13,26)*pmt.stripdx; %(mm) 26-strip LAPPD
    gain=fcnpdfnorm(targets(:,5),targets(:,6),2)./sum(fcnpdfnorm(pmt.stripx,targets(:,6),2),2);  targets(:,4)=targets(:,4).*gain;

%     l=203; %length of anode, 200mm or 800mm;
%     speed = 190; %(mm/ns) pulse speed down anode
%     thr=[.006 .5];  method={'fraction','amplitude'};
%     for i=1
%         t=cell(1,4);  for j=1:4; [~,~,~,tj]=fcnpulsewidth(inputs(:,vi{j})',thr(i),[.006 1E3],[1 1E3],1E6,method{i});  t{j}=tj*.2; end
%         xa = speed*(t{1}-t{2})/2;  ta = (t{1}+t{2} - l/speed)/2; %(mm) and (ns)
%         xb = speed*(t{3}-t{4})/2;  tb = (t{3}+t{4} - l/speed)/2; %(mm) and (ns)
%     end;
%     speed=186;  l=1000;  xc=speed*(ta-tb)/2;  tc=(ta+tb - l/speed)/2; %(mm) and (ns)
%     x=[(xa+xb)/2 xc ta tb tc];  z=targets(:,[1 2 7 8 3]);  e=z-x;  xtof=x;  %fig(2,3,1.5); for i=1:6; sca; plot(x(:,i),z(:,i),'.'); end
%     a=std(e)
    
    Va=zeros(n,512);  v128=1:128;  t=zeros(n,4); tbias=zeros(n,4); %F=lowPassFilter(100E6,200E6,5E9);
    for i=1:4
        I = inputs(:,vi{i});  %I = filter(F,I')';  nb=2^12;  I=round(I*nb)/nb; %12-bit digitized
        [~,~,~,t(:,i)]=fcnpulsewidth(I',.5,[0 1000],[0 1000],1E6,'fraction');  inputs(:,vi{i}) = I;
    end;  t(t==0)=inf;
    for i=1:n   
        ia=round(min(t(i,:)));  ia=floorandceil(ia-10,0,127);  Va(i,:)=inputs(i,ia+[v128 v128+256 v128+256*2 v128+256*3]); tbias(i,3)=ia*.2;
%          ia=round(min(t(i,1)));  ia=floorandceil(ia-55,-256,256);  va=ia+v128;  k=va>0 & va<257; vb=1:sum(k);  %if ia<0; vb=vb-ia; end
%          for j=1:4
%             Va(i,vb+(j-1)*128)=inputs(i,va(k)+(j-1)*256); 
%          end
%          tbias(i,3)=ia*.2;
    end
    Vb=inputs(:,1:2:1024); %PSEC4
    ha=fig(2,1,'19x19cm'); plot(Va(1:10,:)','.-'); sca; plot(Vb(1:10,:)','.-'); fcntight; ha(1).YLim=ha(2).YLim; xyzlabel(ha,'T (5Gs/s sample)','V');
    
%     tic
%     a=zeros(5,4,3); %[layers, gain, type]
%     gv=[.25 .5 1 2]; 
%     sizemethod={'poly1','exp1','rat01'};
%     I=Va;  T=targets(:,[2 3])-tbias(:,[2 3]);
%     for i=1:5
%         hiddenLayers=zeros(1,i);
%         for j=1:4
%             for k=1:3
%                 [i j k]
%                 [~,a(i,j,k)]=NNtrain(I,T,hiddenLayers,3,sizemethod{k},gv(j));
%             end
%         end
%     end
%     toc
    
    
%     I=Va;  T=targets(:,[1 2 3 4])-tbias;  neta=NNtrain(I,T,[0 0 0],1);
%     I=Va;  T=targets(:,1);                netb=NNtrain(I,T,[0 0 0],3);
%     I=Va;  T=targets(:,2);                netf=NNtrain(I,T,[0 0 0],3);
%     I=Va;  T=targets(:,3)-tbias(:,3);     netc=NNtrain(I,T,[0 0 0],3);
%     I=Vb;  T=targets(:,[1 2 3 4]);        netd=NNtrain(I,T,[0 0 0],3);
%     I=Vb;  T=targets(:,4);                nete=NNtrain(I,T,[0 0 0],3);
%     I=[neta(Va')'+tbias, netb(Va')', netc(Va')'+tbias(:,3), netd(Vb')', nete(Vb')', netf(Va')', sum(Vb,2)];  T=targets(:,[1 2 3 4]);   [netg, perf, sigmas]=NNtrain(I,T,[40 20],3); %#ok<ASGLU>
%     s=sprintf('NN.%s',files{fi});  fprintf('Saving ''%s''... ',s);  save('-v6',s,'neta','netb','netc','netd','nete','netf','netg','sigmas','perf');  fprintf('Done.\n\n');
end
%[neth, perf, sigmas]=NNtrain([neta(Va')'+tbias sum(inputs,2)],T,[40 20],1);


j = targets(:,4)<10000;  xa = targets(j,4); T=targets(:,[1 2 3 4]);
x=net(Va')'+tbias;
%x=neth([x sum(inputs,2)]')';
%e=x-T;
%e=netg([neta(Va')'+tbias, netb(Va')', netc(Va')'+tbias(:,3), netd(Vb')', nete(Vb')', netf(Va')', sum(Vb,2)]')' - T;
e = x - T;  std(e)
ha=fig(2,2,'19x19cm');
sca; hist2(xa,e(j,1),100); xyzlabel('PE','X_{anode} (mm)')
sca; hist2(xa,e(j,2),100); xyzlabel('PE','X_{fiber} (mm)')
sca; hist2(xa,e(j,3),100); xyzlabel('PE','T (ns)')
sca; hist2(xa,e(j,4),100); xyzlabel('PE','PE'); ha(4).YLim=[-50 50]; fcntight('xy')

ha=fig(2,2,'19x19cm');
sca; hist2(T(j,1),x(j,1),100); xyzlabel('PE','X_{anode} (mm)')
sca; hist2(T(j,2),x(j,2),100); xyzlabel('PE','X_{fiber} (mm)')
sca; hist2(T(j,3),x(j,3),100); xyzlabel('PE','T (ns)')
sca; hist2(T(j,4),x(j,4),100); xyzlabel('PE','PE'); ha(4).YLim=[-50 50];
% save -v6 FiberNN.NTC7.DRS4.10mV.1mV.0khz.air.20cm.1MeV.mat net
end




function trainNeutronCones
load('/Users/glennjocher/Google Drive/MATLAB/neutrinos/nViewGUI/TSresults/TS FiberGamma NEW4nofitter - 2L SENSL J60035 5V 50pc 512pix EJ-254 150k.mat')
X=MC.xhat;  T=MC.xtrue;
i=find(~all(X==0,2) & ~all(T==0,2) & ~any(isnan(X),2) & ~any(isnan(T),2) & ~any(imag(X)~=0,2) & ~any(T==inf,2) );
X=real(MC.xhat(i,:));  T=MC.xtrue(i,[13 17]); %[angle velocity E0 E0unc]
X=X(:,[1:15 17:end]);

neta=NNtrain(X,T,[40 20],3);  x=neta(X')';  
e=MC.xhat(i,[13 17])-T; e(:,2)=e(:,2)./T(:,2); std(e)
e=x-T; e(:,2)=e(:,2)./T(:,2); std(e)

MC.xhat(i,[13 17])=x; clear X T e x


newMinAE = fcnangle(MC.vertex(i,6:8),MC.xhat(i,6:8)-MC.xhat(i,1:3))*r2d - MC.xhat(i,13);

%fig; fcnhist(cosd(MC.xhat(i,16)),100)
%fig; fcnhist(cosd(newMinAE),100)

MC.xhat(i,16)=newMinAE;  clear X T e x newMinAE
save('/Users/glennjocher/Google Drive/MATLAB/neutrinos/nViewGUI/TSresults/TS FiberGamma NEW4nofitterNN - 2L SENSL J60035 5V 50pc 512pix EJ-254 150k.mat','input','MC','tsv')


end