function results = fcnfastneutron(input,output,handles,PE,G1,plotflag)
MTCflag = ischecked(handles.GUI.realdataflag);
results = [];  nz = sum(output(1).N);  if nz<19 || nz>1E4;  return;  end;  if plotflag; closeallexcept(handles.GUI.figure1); end
input.neutron = load('dEfit.mat');

k = fcnoptimizerk(input,output(1));  k.timeflag=true;  k.reflections=0;  zt = sort(fcnsigmarejection(k.zt(k.ztpid)));
if nanstd(fcnsigmarejection(zt(zt~=0)))>15 || isempty(zt);  return;  end; %timestamps too spread, possible ion feedback event

s = sqrt(k.nhpp);
p0 = s'*k.pxyz/sum(s) * .95;
e0 = fcnguessinitenergy(k,k.nhpp,p0);

mnb = 2; %max number points to try
BIC = inf(mnb,1);
x = cell(mnb,1);
t0 = zt(1);
es = 3; %mm airgap between xhat and detector wall
Lr = input.cube.Lr - es;
A=[];
B=[];
for nb = 1:mnb
    tic
    nx = nb*5;
    
    if nb>1 %WROTE 2013
        A=zeros(nb-1,nx); for i=1:(nb-1); A(i,4+(i-1)*5)=1; A(i,9+(i-1)*5)=-1; end
        B=-ones(nb-1,1)*0.5; %bounces half ns apart
        p0=x{1}(1:3);  t0=x{1}(4);  e0=x{1}(5)/nb;
    end
    x0=zeros(nb,5);  x0(:,1)=p0(1); x0(:,2)=p0(2); x0(:,3)=p0(3); x0(:,4)=t0; x0(:,5)=e0;  x0=reshape(x0',[1 nx]);
    
    LB=[];  UB=[]; Aeq=[]; Beq=[]; tt=sqrt(sum(Lr.^2))*k.ci;
    if input.cube.shapeID~=3 %cylinder=3
        LB=repmat([-Lr zt(1)-tt-5  .01],[1 nb]);  UB=repmat([ Lr zt(end)-tt+30 min(e0*2+.5,10)],[1 nb]);
    end
    %[x{nb}, fx] = fminunc(@(x) fcnfermatpointvectorized(k,x), x0, optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton'));
    %[x{nb}, fx] = fmincon(@(x) fcnfermatpointvectorized(k,x), x0, A, B, Aeq, Beq, LB, UB, [],optimoptions(@fmincon,'Display','off'));
    %[x{nb}, fx] = patternsearch(@(x) fcnfermatpointvectorized(k,x), x0, A, B ,Aeq, Beq, LB, UB, [], optimoptions(@patternsearch,'Display','off','UseVectorized',true,'UseCompletePoll',true));
    %[x{nb}, fx] = simulannealbnd(@(x) fcnfermatpointvectorized(k,x),x0,LB,UB);
    %[x{nb}, fx] = ga(@(x) fcnfermatpointvectorized(k,x),nb*5,A,B,Aeq,Beq,LB,UB,[],optimoptions(@ga,'Display','off','UseVectorized',true));
    %[x{nb}, fx] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
    
    %     if nb==1
    %         [x{nb}, fx] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
    %     else
    %         swarmx=cell(3,1); swarmfx=zeros(3,1);
    %         for i=1:3
    %             [swarmx{i}, swarmfx(i)] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
    %         end
    %         [fx, i] = min(swarmfx);  x{nb}=swarmx{i};
    %     end
    
    input.tsv=5;
    switch input.tsv
        case 1
            return; [x{nb}, fx] = fminunc(@(x) fcnfermatpointvectorized(k,x), x0, optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton'));
        case 2
            [x{nb}, fx] = fmincon(@(x) fcnfermatpointvectorized(k,x), x0, A, B, Aeq, Beq, LB, UB, [],optimoptions(@fmincon,'Display','off'));
        case 3
            [x{nb}, fx] = patternsearch(@(x) fcnfermatpointvectorized(k,x), x0, A, B ,Aeq, Beq, LB, UB, [], optimoptions(@patternsearch,'Display','off','UseVectorized',true,'UseCompletePoll',true));
        case 4
            [x{nb}, fx] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
        case 5
            if nb==1
                [x{nb}, fx] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
            else
                swarmx=cell(4,1); swarmfx=zeros(4,1);
                for i=1:3
                    [swarmx{i}, swarmfx(i)] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true,'SwarmSize',min(100,10*nb*5)));
                end
                [swarmx{4}, swarmfx(4)] = fmincon(@(x) fcnfermatpointvectorized(k,x), x0, A, B, Aeq, Beq, LB, UB, [],optimoptions(@fmincon,'Display','off'));
                [fx, i] = min(swarmfx);  x{nb}=swarmx{i};
            end
    end
    
    BIC(nb) = 2*fx + nx*log(input.cube.pixels+numel(k.ztpid));    %AIC1 = (2*fx + nx*2);    %AICc1 = AIC1 + (2*nx*(nx+1))/(nz-nx-1);
    fprintf('%.0f points (%.2fs): fx=%.1f, BIC=%.1f\n',nb,toc,fx,BIC(nb))
    if nb>1 && BIC(nb)>BIC(nb-1); BIC(nb)=inf; BIC=BIC(1:nb); break; end
end
[~, nb] = min(BIC); fprintf('min BIC found at %.0f points: [',nb); x1 = x{nb}; 
x1 = sortrows(reshape(x1,[5 nb])',4);  %order by time
x1v = reshape(x1', [1 nb*5]); %row vector
fprintf('%5.2f',x1(:,5)); fprintf('] MeV\n');

if plotflag
    [~,~,~,~,g] = fcnfermatpointvectorized(k,x1v);  g=g*sum(x1v(5:5:end));  f=zeros(input.cube.pixels,1);  f(k.pid)=g;
    
    fig(1,2,2);
    sca; popoutsubplot(handles.GUI.axes1,gca);  scatter3(gca,x1(:,1),x1(:,2),x1(:,3),sqrt(x1(:,5))*200,'r','filled'); legend off; fcnview('skew')
    sca; fcnPlotDetector(input,f);              scatter3(gca,x1(:,1),x1(:,2),x1(:,3),sqrt(x1(:,5))*200,'r','filled');
    
%     ha=fig(1,2,1.5,2);  zN=output(1).N;
%     sca; zN(zN==0)=nan; fcnplotdetectorprojection(input,zN>0,zN); title('Measured'); if MTCflag;  title(sprintf('%s event %.0f ',str_(A.pf1),A.events(ei)));  end
%     sca; f(f==0)=nan; fcnplotdetectorprojection(input,f>0,f);  title('Fit')
%     [el,az]=fcnelaz(x1(:,1),x1(:,2),x1(:,3));   [rm,cm]=fcnmapprojection(el*r2d,az*r2d,'winkeltripel');  for i=1:2; plot(ha(i),cm,rm,'b.-','markersize',25); plot(ha(i),cm,rm,'b+','markersize',50,'linewidth',1.5); end
%     fcntight('csigma joint')
end


% if nb==2
%     results.protonprotonflag = [false false];
%     results.minangleerror = 0;
%     results.true = zeros(1,28);
%     results.xhat = results.true;
%     results.true(1:10) = [mean(G1.p1(1:2,:)) mean(G1.t1(1:2,:)) sum(G1.de(1:2)), G1.p1(3,:) G1.t1(3,:) G1.de(3)];
%     results.xhat(1:10) = x1v(1:10);
% end


if nb>1 %at least 3 photons from each bounce
%     %OLD ENERGY METHODS!!
%     dEhat = zeros(2,1);  %pose = zeros(np,8);
%     for i=1:2
%         j=output(1).pid;
%         zN = accumarray(j,sourceProbability(j,i),[input.cube.pixels 1]);
%         
%         %EHAT - MODIFIED Z METHOD MLENERGY
%         dEhat(i) = fcnMLenergy(input,x1(i,1:3),zN);
%         
%         %EHAT - MODIFIED Z METHOD MLPOINT
%         %a = fcnMLpoint(input,zN);  dEhat(i)=a(4); x1(i,1:3)=a(1:3);
%     end
% 
%     %EHAT - MIXTURE POISSON MLPOINT
% %     x0 = zeros(nb,4);  x0(:,4)=exp(nb:-1:1)/exp(nb)*1;  x0(:,1:3) = x1(:,1:3);  x0=reshape(x0',[1 nb*4]);
% %     a = fcnMLpoint(input,output(1).N,x0);  a=reshape(a,[4 nb])';  dEhat=a(1:2,4);  
% %     x1(1:2,1:3) = a(1:2,1:3);
%     
%     %EHAT - MIXTURE POISSON MLENERGY
% %    a = fcnMLenergy(input,x1(:,1:3),output(1).N);  dEhat = a(1:2);
%         
%     x1(1:2,5) = dEhat;

   results = doubleScatterVerification(results,input,G1,PE,handles,x1,ones(size(x1)),'neutron',plotflag);
end

end



function [c,ceq]=nonlcon(x,input)
es = 1; %mm airgap between xhat and detector wall
np = numel(x)/5;
c = zeros(1,np*6);
if input.cube.shapeID==3 %cylinder
    c(1) = norm(x(1:2)) - input.cube.radius+es;
    c(2) = abs(x(3)) - input.cube.height/2+es;
else
    Lr = input.cube.Lr-es;
    for i=1:np
        j = 5*(i-1);
        k = 6*(i-1);
        c(1+k) = abs(x(1+j)) - Lr(1);
        c(2+k) = abs(x(2+j)) - Lr(2);
        c(3+k) = abs(x(3+j)) - Lr(3);
        
        c(4+k) = -x(4+j); %t>0
        c(5+k) = x(4+j)-25; %t<25
        c(6+k) = -x(5+j); %weight>0
    end
end
ceq = [];
j = 5:5:5*np;
ceq(1) = 1 - sum(x(j));
end