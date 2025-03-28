% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [input, output] = fcnBuffer2Memory(input,flags,handles,PE)
MTCflag = ischecked(handles.GUI.realdataflag);

if MTCflag
    [input, D, ZI] = fcnloadMTCQT('',input,flags,handles,1);    %es=%[triggered, amplitude, width, time, integral, prepulse]
    ZI = double(ZI);
    j=find(ZI(:,1));
    if any(j)
        ZF=[j ZI(j,2:end)];
        %load MTCoffsets.mat;  qe=map(:,1).*map(:,2);  ZF(:,2) = ZF(:,2)./(qe(j)+eps)*nanmean(qe); %#ok<NODEF>
                
        gain = 1; %364.83; %nanmean(map(:,2));
        output(1) = oneBookend(ZF,ZI,D,[],0,gain);
    else
        output(1).Nsum=0;
        output(1).D=D;
        output(1).amplitude=0;
    end
else
    gain = input.cube.pmt.amplitude.mu; % (Volts/PE)
    if PE.n==0
        t0 = {0,0};
        pei = {[],[]};
    else
        tspan = input.cube.dsp.t(end);
        t0{1} = 0;                          pei{1} = PE.t>t0{1} & PE.t<tspan;
        t0{2} = (max(PE.t)-tspan*.5);       pei{2} = PE.t>t0{2} & ~pei{1};
    end
    
    for i=1:2
        I = pei{i};  if PE.n>0;  PE.t(I)=PE.t(I)-t0{i};  end
        D = fcnA2D(input,flags,PE,find(I));
        if flags.status.lappd; output(i).Nsum=0; output(i).D=D; continue; end
        %FIXED THRESHOLD NOTE: SENSL J60035 5V 0.100 V amplitude fixed

        [j, a, w, t, integral, Ds] = fcnpulsewidth(D,.5,[input.cube.dsp.noise*6 100],[8 128],input.cube.dsp.dynamicrange(2),'fraction');  t=t*input.cube.dsp.dt+t0{i};
        l = a<input.cube.dsp.noise*6; a(l)=0; w(l)=0; integral(l)=0; t(l)=0;  %these are just noise
        
        ZI = [j a w t integral];
        j=find(a);  ZF=[j ZI(j,2:end)];
        output(i) = oneBookend(ZF,ZI,Ds,I,t0{i},gain); %#ok<*AGROW>
        if PE.n>0 && any(I);  output(i).Ntrue = accumarray(PE.pixel(I), 1, [size(Ds,2) 1]); else output(i).Ntrue=output(i).N*0; end
        
        if i==1 && PE.n>0 && 0
            [pt,jj]=sort(PE.t,'descend');
            T = zeros(input.cube.pixels,1);
            T(PE.pixel(jj)) = pt;
            dt=t-T;  k=t>0 & T<50 & output(1).A>.05;  % & output(1).Ntrue==1;
            DT=dt(k);  A=output(1).A(k);
            hist211(A,DT,100);  xyzlabel('Amplitude (V)','Fit-True Time (ns)')

            
            %NEURAL NETWORK
            k=find(k);  targets=T(k);  inputs=Ds(:,k)';
%             nets=cell(300,1); sigmas=zeros(300,1);
%             for ii=300
%                 net=fcnNNTEST(inputs,targets,[64 16 4]);  nets{ii}=net; %save net.mat net
%                 residuals = targets - net(inputs')';  sigmas(ii)=std(residuals);
%                
%                  load nets93.mat; that=zeros(size(targets));
%                  for ib=1:93
%                      net=nets{ib};  that = that + net(inputs')'/93;
%                      %residuals = targets - net(inputs')';  s(ib)=std(residuals)
%                  end
%             end
            load nets93.mat;  net=nets{71};  residuals = targets - net(inputs')';
            hist211(N,residuals,100);   xyzlabel('Amplitude (V)','NN residual (ns)')
            
            a=[(1:numel(residuals))' residuals];  a=sortrows(a,-2); 
            fig; plot(inputs(a(1:30,1),:)')
            fig; plot(inputs(a(end-30:end,1),:)')

            a=[(1:numel(residuals))' abs(residuals)];  a=sortrows(a,2); 
            fig; plot(inputs(a(1:30,1),:)')
            
            
            sca(handles.GUI.axes1)
        end
        
    end
end



end



function X = fcnA2D(input,flags,PE,i)
pmt = input.cube.pmt;  dsp = input.cube.dsp;
nt = dsp.samples; %number of samples
np = input.cube.pixels; %number of pixels
noise = circshift(input.cube.dsp.N,randi(50),2);

if any(i)
    if flags.status.lappd
        X = zeros(nt, pmt.strips, 2, input.cube.pixels,'single');
        x = [-.5 .5]*input.cube.pixelsize(1);
        [c, strip, pid] = fcnanalogvoltage(input,flags,PE,dsp.t,x,i,true);
        
        c = permute(c,[3 2 1]);
        for j=1:numel(strip)
            X(:,strip(j),1:2,pid(j)) = c(:,:,j);
        end
    else
        X = zeros(nt,np,'single');
        pdf = fcnanalogvoltage(input,flags,PE,dsp.t,0,i,true);
        [pid, ~, J] = fcnunique(PE.pixel(i));
        
        X(:,pid) = fcnaccumrows(pdf,J)';
     end
    X = X + noise;
else
    X = noise;
end
nb = 2^dsp.bits;
X = floorandceil(round(X*nb)/nb,dsp.dynamicrange);

%[V, QT] = fcnphotondeconvolution(input,flags,X);
end


