% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcngettimebases(input,handles,flags)
closeallexcept(handles.GUI.figure1); clc;
[pathname]=uigetdir([pwd '/*.*','MultiSelect','On'],'Select folder:'); 
if isequal(pathname,0); fprintf('No file selected ... Done.\n'); return; end; addpath(pathname); t0=clock;

A = zeros(1536,100,175*6,2,'single'); %T(:,:,:,1)=dt, T(:,:,:,2)=amplitude;
T = zeros(1536,175*6,3);

load('MTCpedestal') %loads X 71x98304
input.MTC.A.pedestals = muvb(8:end,:);

a=dir(pathname);  a=a([a.isdir]);  a=a(3:end);  na=numel(a);
for i=1:na
    [A,T] = process1dir(input,flags,handles,[pathname filesep a(i).name],A,T);
end
save('MTCTimebase.mat','A','T'); fprintf('\nDone. (%gs)',etime(clock,t0))
end


function []=processData()
clc; clear; close all; load('MTCTimebaseS12.mat')  %#ok<*NODEF>
channels = find(any(A(:,:,1,1),2))'; 

X = zeros(3,1536);
for ch = channels
    ch = channels(6)
    x = squeeze(A(ch,:,:,1));  x(x==0)=nan; % fig; plot(x','.')

    
    for sri=1:3; mu=nanmean(x,1); s=nanstd(x,[],1); x(x>mu+2*s | x<mu-2*s)=nan; end;  t=nanmean(x,1); %measured time
    x = squeeze(A(ch,:,:,2));  x(x==0)=nan;  for sri=1:3; mu=nanmean(x,1); s=nanstd(x,[],1); x(x>mu+2*s | x<mu-2*s)=nan; end;  a=nanmean(x,1); %measured amplitude
    CI = fcncol(ones(175,1)*(1:6));%Coarse Index 0-5
    
    %x=mod(x+64,128);
    di=[false abs(diff(t))<2];  t=t(di); %cut jumps
    z = T(ch,:,2)'/1000 + T(ch,:,3)'*7.875; 
    z=z(di); CI=CI(di); 
    [~,i]=min(diff(t));  GI=z'>z(i) & t<t(i); % plot(x(GI),y(GI),'.'); %Gap Index 0-1
    
    %LLS  x = [m b CoarseGap]
    H = zeros(numel(t),3);  H(:,1)=t+128*GI;  H(:,2)=1; for i=1:6; H(CI==i,3)=i; end;
    X(:,ch) = (H'*H)\H'*z;
    residuals = z-H*X(:,ch);  %residuals = residuals - residuals(find(GI,1,'first'));
    
    ha=fig(3,1,'29x19cm'); %plot(x,y,'.');
    for i=1:6
        j=CI==i;
        plot(ha(1),t(j),z(j),'.','Display',['coarse ' num2str(i)]); 
        plot(ha(2),t(j),residuals(j),'.','Display','residuals');
        plot(ha(3),t(j),a(j),'.','Display','amplitude');
    end
    plot(ha(1),t,H*X(:,ch),'b-','Display','fit'); legend('Location','Best')
    [S, R, C] = MTCpixelID2SRCCH(ch);
    xyzlabel(ha,'T (bin)','T (ns)'); ylabel(ha(2),'residuals (ns)'); ylabel(ha(3),'amplitude (bin)'); ha(1).Title.String=sprintf('S%g R%g C%g P%g - %.0fps/bin',S,R,C,ch,X(1,ch)*1000);
    fcnmarkersize(15); fcnfontsize(16); fcntight; ha(3).YLim(1)=0;  for i=1:numel(ha); ha(i).YLabel.Position(1)=-9; end
end


end


function [A,T] = process1dir(input,flags,handles,pathname,A,T)
a = dir(pathname); 
a=a(~cellfun(@isempty,strfind({a(:).name},'.glenn')));
%pattern = 'timebase_scan_fineDelay_(\d*)_coarseDelay_(\d*)_rowcolch_(\d)(\d)(\d)_';
pattern = 'timebasescan_delaynum_(\d*)_fineDelayActual_(\d*)_coarseDelay_(\d*)_rowcolch';

for i=1:numel(a);
    f = a(i).name;
    t = regexpi(f, pattern, 'tokens'); t=t{1};  fine=eval(t{2}); coarse=eval(t{3}); ti=eval(t{1})/5+coarse*175+1;
    
    input = fcnloadMTCQT([pathname filesep f],input,flags,handles,0);  if isempty(input.MTC.C); continue; end
    x=input.MTC.C(:,2);  nx=numel(x);  y=cell2mat(x); %event data
    yi=reshape(y(:,4),[1536 nx]);  j=any(yi,2);  A(j,:,ti,1) = full(yi(j,:)); %dt
    yi=reshape(y(:,2),[1536 nx]);  j=any(yi,2);  A(j,:,ti,2) = full(yi(j,:)); %amplitude
    T(j,ti,:) = [ti fine coarse]; %time index, fine delay, coarse delay
end
end


function plotresiduals()
slopes = zeros(1536,5);  stdres=slopes;  meanamplitudes=slopes;
delay = [0 8 16 24 32];  A=cell(1536*5,2);
for i=1:5
    load(sprintf('CableDelay%.0f',delay(i)))
    for pixel = 1:1536
        pixel
        U = squeeze(T(pixel,:,:,1));
        x = repmat(offsets,[size(U,1) 1])/100;  j=U~=0;  x=x(j); U=U(j);
        
        if i==3
            U(U>100) = U(U>100)-128;
        end

        if any(U)
            ft=fittype('poly1');  fr=fit(x, U, ft);  slopes(pixel,i)=1/fr.p1;
            %hist211(x,U); xyzlabel(sprintf('Pixel %.0f offset (ns)',pixel),'T (128-vector bin)');
            %xa = linspace(min(x),max(x),100); plot(xa,fr(xa),'k','linewidth',2);
            residuals = U - fr(x);  stdres(pixel,i)=std(residuals);  A{pixel+(i-1)*1536,1}=residuals;
            amplitudes = T(pixel,:,:,2);  amplitudes=amplitudes(j);  A{pixel+(i-1)*1536,2}=amplitudes;
            meanamplitudes(pixel,i)=mean(amplitudes);
            
            %hist211(x,residuals); xyzlabel(sprintf('Pixel %.0f offset (ns)',pixel),'T (128-vector bin)');
        end
    end
end
B=cell2mat(A);  i=abs(B(:,1))<2 & B(:,2)<1300 & B(:,2)>100;  hist211(B(i,1),B(i,2)); xyzlabel('t residual (bins)','pulse height (bins)')
i = stdres~=0;  hist211(stdres(i),meanamplitudes(i)); xyzlabel('residual std (bins)','mean pulse height (bins)')



fxva = mean(meanamplitudes,2);
fig; fcnplotdetectorprojection(input,fxva~=0,fxva,[0 0 0],'winkeltripel','SCROD')
h=colorbar; h.Label.String = 'pulse heights (bins)';

load('MTCmapping.mat');
[double(X(:,4)), slopes(X(:,1),:)]

end


% function combinefiles()
% %COMBINE ALL FILES --------------------------------------------------------
% delay = [0 8 16 24 32];
% clear U Ua x xa
% pixel = 260;
% 
% Ua = []; xa=[];
% for i=1:5
%     load(sprintf('CableDelay%.0f',delay(i)))
%     
%     U = squeeze(T(pixel,:,:));
%     x = repmat(offsets,[size(U,1) 1])/100;
%     
% %     if i>2; 
% %         x=x-15;
% %     end
%     
%     Ua = [Ua; U];
%     xa = [xa; x+delay(i)];
% end
% 
% x=xa; U=Ua; i=U~=0;
% %U = U-x/.3666;
% if any(i);  hist211(x(i),U(i)); xyzlabel(sprintf('Pixel %.0f offset (ns)',pixel),'T (128-vector bin)'); end
% end
