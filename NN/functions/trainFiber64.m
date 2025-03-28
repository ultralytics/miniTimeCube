% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function trainFiber64(files)
if nargin==0
    applyFlag=false;
    %files=uigetfile('*.mat','MultiSelect', 'on');  if ischar(files); files={files}; end
    files={'TS.25ns FiberPointSource - 1L SENSL J60035 5V 128ch BCF10 100k.mat'};
else
    applyFlag=0;
end

close all
nf=numel(files);
for fi=1:nf
    fprintf('Loading ''%s''... ',files{fi});  load(files{fi});  fprintf('Done.\n');
    inputs=MC.xhat;  T=MC.xtrue;  input.extra=[]; clear tsv MC %T=targets, I=inputs
    
    vi{1}=1:256; vi{2}=vi{1}+256;  a1=max(inputs(:,vi{1}),[],2);  a2=max(inputs(:,vi{2}),[],2);
    goodAmplitude=all([a1 a2]>(input.cube.pmt.amplitude.mu+input.cube.dsp.noise*0)*10,2);
    goodRatio=max(a1./a2,a2./a1)<10;
    i=find(~all(inputs==0,2) & ~all(T==0,2) & ~any(isnan(inputs),2) & goodAmplitude & goodRatio & T(:,3)>10);  i=i(1:end); inputs=inputs(i,:);  T=T(i,:);  n=size(inputs,1);  numel(i)
    
    t=zeros(n,2);  a=t; w=t; %F=lowPassFilter(100E6,200E6,5E9);
    for i=1:2
        %I = inputs(:,vi{i}); I = filter(F,I')';  nb=2^12;  I=round(I*nb)/nb;  inputs(:,vi{i})=I;  %12-bit digitized
        [~,a(:,i),w(:,i),t(:,i),~]=fcnpulsewidth(inputs(:,vi{i})',.5,[0 1000],[0 1000],1E6,'fraction');
    end
    ci=abs(diff(t,1,2))<50; inputs=inputs(ci,:); T=T(ci,:); I=inputs;  
    In = I./std(I,[],2) - mean(I,2);  % normalize rows
    fig; plot(I(1:8,:)','.-')
    
    name=sprintf('NN.%s',files{fi});  % good ratio %improved ratio %2normalizations
    % outputs = T; save wavedata3ns.mat inputs outputs name
    if applyFlag
        load(name);
        pos=net{1}(In')';  Id=[NNwaveformStats(I(:,1:256),2) NNwaveformStats(I(:,257:end),2) pos(:,1)];
    else
        [net{1}, perf(1), sigmas(1:2)] = Copy_of_NNtrain(In,T(:,[1 2]),32,1);
        pos=net{1}(In')';  Id=[NNwaveformStats(I(:,1:256),2) NNwaveformStats(I(:,257:end),2) pos(:,1)];   [net{2}, perf(2), sigmas(3)]  =  NNtrain(Id,T(:,4),[20 20],1);
    end
    xhat(1:2) = conventionalFiberFit(inputs, T, [.006, 0.5]);
    xhat{3} = [pos net{2}(Id')'];
    
    %PLOT RESULTS
    hf=findobj(0,'Name','NNplot');  if isempty(hf);  [ha,hf]=fig(2,3,'20x30cm');  hf.Name='NNplot';  else;  figure(hf);  end
    nb=100; clear S;  sx={'error (mm)','error (ns)','E (MeV)'};  st={'Position','Time','Energy'};  sd={'CA 6 mV','CF 0.5','WAVE'}; ha(2).Title.String=name;
    ci=T(:,4)>0 & T(:,4)<1.0;   
    xa=T(ci,4);  xlabel = 'E (MeV)';
    for i=1:numel(xhat)
        e = xhat{i} - T(:,[1 2 4]);
        nanstd(e)
        for j=1:3
            sca;
            if(i)<3 && j==3; continue; end
            if(i)==1; cla; end
            y=e(ci,j);  hist2(xa,y,nb,'scatter');   [~,~,S{i,j}] = movingMean(xa,y,nb,0);  xyzlabel(xlabel,sx{j},'',st{j});
        end
        for j=1:3
            sca;
            if(i)<3 && j==3; continue; end
            x=linspace(min(xa),max(xa),nb); plot(x,S{i,j}.s(x),'.-','Display',sd{i});  xyzlabel(xlabel,sx{j});
        end
    end
    fcntight(ha(4:6),'y0'); fcnmarkersize(ha(1:3),5); fcnfontsize(16);
    
    %% TESTING
    [~,i1]=max(inputs(:,vi{1}),[],2);  [~,i2]=max(inputs(:,vi{2}),[],2);
    fig; hist2(xhat{3}(:,2),i2,30)
    
    i = round((xhat{3}(:,2)/0.200)); Ic=inputs;
    for j=1:size(inputs,1)
        %Ic(j,:) = [circshift(inputs(j,1:255),-i(j)), circshift(inputs(j,256:end),-i(j))];
        Ic(j,:) = circshift(inputs(j,:),-i(j));
    end
    Ic = Ic./std(Ic,[],2) - mean(Ic,2);  % normalize rows
    fig; plot(Ic(1:8,:)','.-')
    %xhat{4}=[net{1}(Ic')'+[i*0 i*.2 i*0] net{2}(Id')'];
    that = net{1}(Ic')'; that=that(:,2)+i*.2;
    std(that-T(:,2))
    %% TESTING DONE
    
    
%     %HYPERPARAMETER STUDY (amplitude and fraction)
%     n = 500;
%     fractions = linspace(0.07,.9999,n);
%     amplitudes = linspace(0.004,.2,n);
%     truth = T(:,[1 2 4]);
%     y = zeros(n,12);
%     for i=1:n
%         i
%         a = conventionalFiberFit(inputs, T,[amplitudes(i), fractions(i)]);
%         b = nanstd(a{1}-truth); eb = 1-nanmean(isnan(a{1}),1);
%         c = nanstd(a{2}-truth); ec = 1-nanmean(isnan(a{2}),1);
%         y(i,:) = [b, eb, c, ec];
%     end
%     fig(2,1,'19cm','19cm',[2, 2, 1.9, 1.9, 1.9, 1.2])
%     [h,h1,h2] = plotyy(amplitudes,y(:,1),amplitudes,y(:,2));  fcnlinewidth(3);
%     title('Constant Amplitude (CA) hyperparameter study'); xlabel('threshold (mV)'); ylabel(h(1),'error 1\sigma (mm)'); ylabel(h(2),'error 1\sigma (ns)'); fcnfontsize(16); fcntight('x0')
%     c=[0 0.447 0.741]*.6;  h1.Color=c;  h(1).YColor=c; c=fcnlightencolor([0.929  0.694  0.125],.3);  h2.Color=c;  h(2).YColor=c;
%     sca; plot(amplitudes,y(:,4),'.-'); xyzlabel('threshold (mV)','candidates (fraction)')
%     
%     fig(2,1,'19cm','19cm',[2, 2, 1.9, 1.9, 1.9, 1.2])
%     sca; [h,h1,h2] = plotyy(fractions,y(:,7),fractions,y(:,8));  fcnlinewidth(3);
%     title('Constant Fraction (CF) hyperparameter study'); xlabel('threshold (fraction)'); ylabel(h(1),'error 1\sigma (mm)'); ylabel(h(2),'error 1\sigma (ns)'); fcnfontsize(16); fcntight('x0')
%     c=[0 0.447 0.741]*.6;  h1.Color=c;  h(1).YColor=c; c=fcnlightencolor([0.929  0.694  0.125],.3);  h2.Color=c;  h(2).YColor=c;
%     ha=sca; plot(fractions,y(:,10),'.-'); xyzlabel('threshold (fraction)','candidates (fraction)'); ha.YLim=[0 1.01];
    
    %SAVE
    if ~applyFlag;  fprintf('Saving ''%s''... ',name);  save(name,'net','sigmas','perf','name','S');  fprintf('Done.\n\n');  end
end
