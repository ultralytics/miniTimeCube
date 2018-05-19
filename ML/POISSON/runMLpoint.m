function results = runMLpoint(input,output,handles,G1,plotflag)
results=[];  MTCflag = ischecked(handles.GUI.realdataflag); if MTCflag; A=input.MTC.A; ei=input.eventNumber+1; end; if plotflag; try closeallexcept(handles.GUI.figure1); catch; end; end; tic

%EVENT REJECTION
nz = output(1).Nsum;  zt = output(1).t;
if ~any(zt) || nanstd(fcnsigmarejection(zt(zt~=0)))>15 || nz<3 || nz>1E4;  return;  end; %timestamps too spread, possible ion feedback event

%FITTERS [XYZ T E]
type = 'Poisson';
switch type
    case 'Poisson'
        X = [fcnMLpoint(input,output(1)) 0]; %[x y z e]
        results.xhat = [X(1:3) 0 X(4)]; %[xyzte]
    case 'Time'
        x=cell(2,1);
        for i=1:1
            t = output(i).t;  t=t(t~=0);
            if numel(t)>4
                k = fcnoptimizerk(input,output(i));  k.reflections=0;  k.timeflag=1;
                
                s = sqrt(k.nhpp);
                p0 = s'*k.pxyz/sum(s) * .95 * 0;
                t0=min(fcnsigmarejection(t,2,3));
                e0 = fcnguessinitenergy(k,k.nhpp,p0);
                
                nb = 1;  nx = nb*5;
                x0 = zeros(nb,5);  x0(:,5)=e0/nb;  x0(:,4)=t0;  x0(:,1:3)=p0;  x0=reshape(x0',[1 nx]);
                es = 1; %mm airgap between xhat and detector wall
                Lr = input.cube.Lr - es;
                
                A=[];  B=[];  LB=[];  UB=[];  Aeq=[];  Beq=[];
                if input.cube.shapeID~=3
                    LB=repmat([-Lr t0-sqrt(sum(Lr.^2))*k.ci-10  .001],[1 nb]);  UB=repmat([ Lr t0+30 10],[1 nb]);
                end
                
                %[x{i}, fx] = fmincon(@(x) fcnfermatpointvectorized(k,x), x0, A, B, Aeq, Beq, LB, UB, [],optimoptions(@fmincon,'Display','off'));
                %[x{i}, fx] = patternsearch(@(x) fcnfermatpointvectorized(k,x), x0, A, B ,Aeq, Beq, LB, UB, [],optimoptions(@patternsearch,'Display','off','UseVectorized',true,'UseCompletePoll',true));
                %[x{i}, fx] = simulannealbnd(@(x) fcnfermatpointvectorized(k,x),x0,LB,UB);
                %[x{i}, fx] = ga(@(x) fcnfermatpointvectorized(k,x),nb*5,A,B,Aeq,Beq,LB,UB,[],optimoptions(@ga,'Display','off','UseVectorized',true));
                [x{i}, fx] = particleswarm(@(x) fcnfermatpointvectorized(k,x),nb*5,LB,UB,optimoptions(@particleswarm,'Display','off','UseVectorized',true));
            end
        end
        X=x{1};
        results.xhat=X;
end
toc

%TRUTH
w=G1.de.*double(G1.p1inside);
results.true = [weightedMean((G1.p1+G1.p2)/2,w) 0 sum(w)];  

if plotflag
    k=fcnoptimizerk(input,output(1)); Xf=[X(1:3) 0 1];  k.reflections=0;  %[xyztw];
    [~,~,~,~,g] = fcnfermatpointvectorized(k,Xf);  g=g*results.xhat(5);  f=zeros(input.cube.pixels,1);  f(k.pid)=g;
    
    ha=fig(1,2,1.5,2);  zN=output(1).N;
    sca; zN(zN==0)=nan; fcnplotdetectorprojection(input,zN>0,zN); title('Measured'); %if MTCflag;  title(sprintf('%s event %.0f ',str_(A.pf1),A.events(ei)));  end
    sca; f(f==0)=nan; fcnplotdetectorprojection(input,f>0,f);  title('Fit')
    [el,az]=fcnelaz(X(1),X(2),X(3));   [rm,cm]=fcnmapprojection(el*r2d,az*r2d,'winkeltripel');  for i=1:2; plot(ha(i),cm,rm,'b.-','markersize',25); plot(ha(i),cm,rm,'b+','markersize',50,'linewidth',1.5); end
    %sca; fcnPlotDetector(input,f); plot3(X(1),X(2),X(3),'b.-','markersize',25); title('Fit')
    fcntight('csigma joint')
end