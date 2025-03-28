% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function results = fcnmuon(input,output,handles,photons,G1,plotflag)
results=[];  MTCflag=ischecked(handles.GUI.realdataflag);  if plotflag; closeallexcept([handles.GUI.figure1, findobj('Name','Event Viewer')]); end
nz = output(1).Nsum;  if nz<30 || nz>30E4;  return;  end
k = fcnoptimizerkr(input,output(1));  k.reflections=0;  nppl = 20; %number of points per line
zN = output(1).N;

zt = output(1).t;  minzt = min(zt)*0;
k.timeflag = false;


%INITIAL GUESS
igmethod = 'center'; %initial guess method
switch igmethod
    case 'brightestpixel'
        [~,i] = max(zN);  p1 = input.cube.all.xyz(i,:);
    case 'firstlight'
        xyzt = sortrows([input.cube.all.xyz(output(1).pid,:) output(1).t],4);
        p1 = mean(xyzt(1:min(nz,50),1:3));
    case 'center'
        p1 = fcnMLpoint(input,zN);  p1=p1(1:3);
end
[~, az, el] = fcnuniformsphere(20,40);  np = size(az,1);  t0 = repmat(minzt,[np 1]);  p0 = repmat(p1,[np 1]);

k.Xrtype = 2;  
X = reduced2full(input,k,[p0 t0 el az]);  fx0 = fcnfermatpointvectorized(k,zt,X);  [~,i] = min(fx0);
k.Xrtype = 2;% 1-[xyz xyz t] 7DOF,    2-[xyz t elaz] 6DOF,    3-[elaz elaz t] 5DOF
X0r = full2reduced(k,X(i,:));

%OPTIMIZE
Lr = input.cube.Lr - 3;  %mm airgap between xhat and detector wall
A=[];  B=[];  LB=[];  UB=[];  Aeq=[];  Beq=[];  tlb = minzt-sqrt(sum(Lr.^2))*k.ci;  tub = minzt+30;
if input.cube.shapeID~=3  %cylinder=3
    switch k.Xrtype
        case 1 %[xyz xyz t]
            LB=[-Lr -Lr tlb];               UB=[Lr Lr tub];
        case 2 %[xyz t elaz]
            LB=[-Lr tlb  -pi/2 -pi];        UB=[Lr tub pi/2 pi];
        case 3 %[elaz elaz t]
            LB=[-pi/2 -pi -pi/2 -pi tlb];   UB=[pi/2 pi pi/2 pi tub];
    end
end
[Xr, fx] = patternsearch(@(x) fermatLineVectorized(input,k,zt,x), X0r, A, B ,Aeq, Beq, LB, UB, input.optimizer.psoptions3);
%[Xr, fx] = fmincon(@(x) fermatLineVectorized(input,k,zt,x), X0r, A, B ,Aeq, Beq, LB, UB,[],input.optimizer.options2);
%BIC = 2*fx + nx*log(nz);    %AIC1 = (2*fx + nx*2);    %AICc1 = AIC1 + (2*nx*(nx+1))/(nz-nx-1);

k.timeflag = true;
switch k.Xrtype
    case 1 %[xyz xyz t]
        Xr2 = Xr([4:6, 1:3, 7]);
    case 2 %[xyz t elaz]
        Xr2 = [Xr(1:4) fcnelaz(-fcnelaz2CC(Xr(5:6)))];
    case 3 %[elaz elaz t]
        Xr2 = Xr([3:4, 1:2, 5]);
end
results.fx1 = fermatLineVectorized(input,k,zt,Xr);  %try forward
results.fx2 = fermatLineVectorized(input,k,zt,Xr2); %try backward
if results.fx2<results.fx1;  Xr=Xr2;  end %flip direction!


X = reduced2full(input,k,Xr);
xreshaped = reshape(X,[5 nppl])';
pos = xreshaped(:,1:3);

%ENERGY
rhat = rangec(pos(1,:)/.95,pos(end,:)/.95); %mm
dEhat = mean(fcnguessinitenergy(k,zN,pos)); %MeV
dE = sum(G1.de(G1.p1inside & G1.parentinside)); %MeV

%RESULTS
if MTCflag
    results.true = [0 0 0 0];
else
    i = find(G1.pid==13,1,'first');  
    if any(i)
        p1=G1.p1(i,:);  vec = G1.p2(i,:)-p1;    elaz = fcnelaz(vec); %rad
        pn = fcndetectorintercept(input,p1,vec);  r = rangec(p1,pn);
        results.true = [elaz, dE, r];
        results.angleerr = acosd( dot((pn-p1)/r, fcnvec2uvec(pos(end,:)-pos(1,:))) );
    end
end
elazhat = Xr(5:6);
results.xhat = [elazhat, dEhat, rhat];

if plotflag
    xi=1:5:nppl*5;  yi=xi+1;  zi=xi+2;  k = fcnoptimizerk(input);  k.timeflag = false;
    [~, ~, fall] = fcnfermatpointvectorized(k,zt,X); fall = fall./k.QEmap;
    a=sum(fall,2); a=a/max(a); b=a; %zeros(1536,1);  b(k.upid)=a;
    
    ha=fig(1,2,2);
    sca;  zN(zN==0)=nan;  fcnplotdetectorprojection(input,zN>0,zN); title('Measured  z')
    sca;  b(b==0)=nan;  fcnplotdetectorprojection(input,b>0,b);  title('Reconstructed  z=H(x)')
    [el,az]=fcnelaz(X(xi),X(yi),X(zi));   [rm,cm]=fcnmapprojection(el*r2d,az*r2d,'winkeltripel');  for i=1:2; plot(ha(i),cm,rm,'b.-','markersize',25); plot(ha(i),cm(end),rm(end),'b+','markersize',50,'linewidth',1.5); end
    
%     popoutsubplot(handles.GUI.axes1, fig(1,1,2))
%     h=plot3(X(xi),X(yi),X(zi),'b.-'); h.MarkerSize = 30;
end
end


function fx = fermatLineVectorized(input,k,zt,xr)
x = reduced2full(input,k,xr);
fx = fcnfermatpointvectorized(k,zt,x);
end


function X = reduced2full(input,k,Xr)
switch k.Xrtype
    case 1 %[xyz xyz t]
        p1 = Xr(:,1:3);
        p2 = Xr(:,4:6);  T = Xr(:,7);
    case 2 %[xyz t elaz]
        vel=fcnelaz2CC(Xr(:,5:6));
        p1 = fcndetectorintercept(input,Xr(:,1:3),-vel);
        p2 = fcndetectorintercept(input,Xr(:,1:3), vel);  T = Xr(:,4);
    case 3 %[elaz elaz t]
        p1 = fcndetectorintercept(input,0,fcnelaz2CC(Xr(:,1:2)));
        p2 = fcndetectorintercept(input,0,fcnelaz2CC(Xr(:,3:4)));  T = Xr(:,5);
end
c = 299.792458; %speed of light (mm/ns);
nppl = 20;
dpv = linspace(0,1,nppl); %delta points vector
[r,dx] = rangec(p2,p1); %r (mm)
X = zeros(size(Xr,1),5*nppl);
for i=1:3
    X(:,i:5:5*nppl) = p1(:,i)-dx(:,i)*dpv;
end
X = X*.95; %keep away from walls
X(:,5:5:nppl*5) = 1/nppl; %weighti
X(:,4:5:nppl*5) = T + (r/c)*dpv; %t (ns)
end


function Xr = full2reduced(k,X) %ONLY WORKS FOR 1, NOT VECTORIZED!!!!
[~,np]  = size(X);  np = np/5; %nv = number of vectorized points
posi    = (1:5:np*5)' + [0 1 2];
timei   = 4:5:np*5;
pos     = X(posi);  p1=pos(1,:);  p2=pos(end,:);
t       = X(timei(1));
switch k.Xrtype
    case 1 %[xyz xyz t]
        Xr = [p1 p2 t];
    case 2 %[xyz t elaz]
        Xr = [(p1+p2)/2 t fcnelaz(p2-p1)];
    case 3 %[elaz elaz t]
        Xr = [fcnelaz(p1) fcnelaz(pos2) t];
end
end


