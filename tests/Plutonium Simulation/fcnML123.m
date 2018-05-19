function [newHandles, MC, d] = fcnML123(input,handles,flags,d,table)
newHandles = []; MC = [];  %flags.status.CRLB=1;
if input.dValid.count<1 || (flags.status.ML1==0 && flags.status.ML2==0 && flags.status.ML3==0)
    deleteh(handles.ML1);
    deleteh(handles.ML2);
    deleteh(handles.ML3);
    return
end
flags.status.CRLB=1;
fprintf('Running ML123...'); startclock = clock; %#ok<NASGU>
om1     = ones(input.nxy^2, input.nrp);
L       = zeros(input.nxy^2, input.nrp);
ov1     = ones(input.nxy^2, 1);
nbatch  = 100; %number of measurments per batch max
sk      = table.mev.de/4/pi; %stabilizing constant
pdfur   = zeros(input.nxy^2, nbatch);
vnrp    = 1:input.nrp;
e       = input.fluxnoise.systematic.estimated;

nv      = input.dValid.count;
Civ     = zeros(nv,1);
Tiv     = zeros([size(L) nv]);
rpscale = input.rpVec;  irpscale=1./rpscale;
for iv = 1:nv
    id=input.dValid.idx(iv);  d1=d(id);  fprintf(' D%.0f...',d1.number)
    
    [Bi, ~, d1]=fcnenergycut(input, flags, table, d1);  urtable_epdf = d1.est.urtable.epdf*sk; %maintain numerical integrity
    Si = interp1(d1.est.urtable.r, d1.est.urtable.n, d1.est.r) * rpscale; %number of events from unknown source
    Ti = Bi+Si; %Bi = # background events, Ti = # total events
    Ci = input.dValid.validevents(iv); %Ci = # candidate measurements
    if flags.status.CRLB && ~flags.status.ML3
        Ci = d1.n.all;
    end

    %ML3-------------------------------------------------------------------
    if flags.status.ML3 && Ci>0
        Ci = input.dValid.validevents(iv); %Ci = # candidate measurements
        ae1=fcnprob3(input, table, d1, flags);  d1.est.ae1=ae1;  
        allb = (ae1.kr*e(1) + ae1.mantle*e(2) + ae1.crust*e(3) + ae1.fn*e(4) + ae1.acc*e(5) + ae1.cosm*e(6))*sk;
        
        nb = ceil(Ci/nbatch); %number of batches
        i2 = fcnindex1c(d1.est.urtable.r, d1.est.r, '*exact');  i1=floor(i2);  i3=ceil(i2);  f=i2-i1; %range interpolant coefficients
        for ib  = 1:nb %number of 100-measurement batches
            v0 = (nbatch*ib-nbatch+1):min(nbatch*ib,Ci);  v1=1:numel(v0);  ei=d1.z.ei(v0);  Lb=allb(v0);  Lc=om1;
            ctip = fcnindex1c(d1.snr.ct, d1.est.puvec*d1.z.udxecef(v0,:)'); %cosine theta params
            for i = v1
                ev=urtable_epdf(:,ei(i));  e1=ev(i1);  epdf=e1+f.*(ev(i3)-e1);
                pdfur(:,i) = epdf.*d1.snr.apdf(ctip(:,i));
            end
            for j = vnrp
                Lbj = Lb*irpscale(j);
                Lc1=ov1;
                for i = v1
                    Lc1 = Lc1.*(Lbj(i)+pdfur(:,i));
                end
                Lc(:,j)=Lc1;
            end
            x=fcnminmax(Lc);
            L = L + log(Lc/sk); %fcnminmax(Lc)
        end
        g = ov1*irpscale; %normalize denominator
        L = L - log(Ti.*g)*Ci;
    end
    
    
    %ML2-------------------------------------------------------------------
    if flags.status.ML2 && Ci>0
        zn = accumarray(d1.z.eic, 1, [table.mev.ne 1]); %measurments per bin;
        epdf = e(1)*d1.epdf.kr + e(2)*d1.epdf.mantle + e(3)*d1.epdf.crust + e(4)*d1.epdf.fastneutron + e(5)*d1.epdf.accidental + e(6)*d1.epdf.cosmogenic;
        if flags.status.CRLB
           zn = d1.epdf.all*table.mev.de;
           epdf = d1.epdf.allbackground;
        end
        
        zei=find(zn);  zn=zn(zei);  ne=numel(zei);
        nr = numel(d1.est.urtable.n);
        bpdf = ones(nr,1)*epdf(zei); %background pdf
        urpdf = d1.est.urtable.epdf(:,zei);
        
        a = zeros(nr, input.nrp, ne);
        for i = 1:input.nrp
            a(:,i,:) = urpdf*rpscale(i) + bpdf;
        end
        a = log(a);
        
        b = zeros(nr, input.nrp);
        for i = 1:ne
            b = b + a(:,:,i)*zn(i);
        end
        
        L = L + interp1(d1.est.urtable.r, b, d1.est.r) - log(Ti)*Ci;
    end

    
    %ML1-------------------------------------------------------------------
    Tiv(:,:,iv)=Ti;  Civ(iv)=Ci;
    
    
    d(id)=d1;
    %D1--------------------------------------------------------------------
    if iv==1
        L1=L;
    end
end

C = fcndcorr(input, d); %cov mat
Ln = lognormal(Civ(1), Tiv(:,:,1), sqrt(C(1,1))); L1 = L1+Ln;
if nv>1
   [~,Ln]=fcnmultinormpdf(Tiv,Civ,C);  L=L+Ln;
else
   L = L1; 
end
%L = fcnrpprior(input, flags, L); %reactor power prior

MLplot
end


function p = lognormal(x,u,s)
p = - log(sqrt(2*pi)*s) - (x-u).^2./(2*s.^2);
end


function [C, D] = fcndcorr(input, d)
%Count Covariance from one detector to another
%C = covariance matrix
%D = correlation matrix

C = ones(input.dValid.count);
for i = 1:input.dValid.count
    d1 = d(input.dValid.idx(i));   

    Bi = d1.n.allbackground; %mean background rate
    C(i,:) = C(i,:)*d1.est.br1s;
    C(:,i) = C(:,i)*d1.est.br1s;
    C(i,i) = C(i,i) + d1.est.d1s.^2 + Bi + Bi*.025; %Bi*.025 adds osc
end

if nargout>1
    D = ones(input.dValid.count);
    for i = 1:input.dValid.count
        D(:,i) = C(:,i)./(C(i,i)*diag(C)).^.5;
    end
end

end
