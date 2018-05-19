function [fall, r] = fcnsolidanglevectorized(k,X)
%mode = 1; %1=simple, 2=reflect
np = size(X,1);

switch k.reflections
    case 0 %do nothing
        [fall,r] = solidanglec([X zeros(np,3)],k.pxyz,k.nxyz,k.pa,k.al,k.yFt); %(pixelxyz(mm),pixelnormalxyz(uvec),pointxyz(mm),pixelarea(mm^2),attenuationlength(mm))
    otherwise %reflection model
        F = cell(np,1);  R=F;  Lr=k.Lr*2;
        for i=1:np
            b=[X(i,1:3) 0 0 0 0];  a=b; %[x y z nxr nyr nzr gen]
            for g=1:k.reflections %gen
                a = fcnreflectedpoints(Lr,a);  a(:,7)=g;  b = [b; a];
                %[~,j] = fcnunique(round(b(:,1:3)*1E3)/1E3,'rows');  b = b(sort(j),:);  a = b(b(:,7)==g,:);
                [~,j] = fcnunique(round(b(:,1:3)*1E3)/1E3,'rows');  b = b(sort(j),:);  a = b(b(:,7)==g,:);
                %if g==1; h=[]; popoutsubplot(evalin('base','handles.GUI.axes1'),fig); end; h(g)=fcnplot3(a,'.','markersize',25/sqrt(g),'color',fcndefaultcolors(g),'DisplayName',sprintf('%.0f gen%.0f pts',size(a,1),g)); axis tight; legend(h)
            end
            [F{i}, R{i}] = solidanglec(b,k.pxyz,k.nxyz,k.pa,k.al,k.yFt);
        end
        s = [k.npixels np size(b,1)];
        fall = reshape(cat(1,F{:}),s);        
        r = reshape(cat(1,R{:}),s); 
end

fall = fall.*k.QEmap;

end



function b = fcnreflectedpoints(detectorWidth,P)
ov = ones(6,1);
np = size(P,1);

%xrv = [1 0 0 1 0 0];
%yrv = [0 1 0 0 1 0];
%zrv = [0 0 1 0 0 1];

xyzrv = [1 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 0 0 1]';

b = [];
for i = 1:np
    p = P(i,1:3);  %[x y z face gen parentid id]
    a = p(ov,:);
    rp =  detectorWidth - p; %reflected positive
    rn = -detectorWidth - p; %reflected negative
    
    a(1,1) = rp(1); %+x  1. FRONT
    a(2,2) = rp(2); %+y  2. RIGHT
    a(3,3) = rp(3); %+z  6. DOWN
    a(4,1) = rn(1); %-x  3. BACK
    a(5,2) = rn(2); %-y  4. LEFT
    a(6,3) = rn(3); %-z  5. TOP
    
    %a(:,4) = P(i,4) + xrv; %x reflections
    %a(:,5) = P(i,5) + yrv; %y reflections
    %a(:,6) = P(i,6) + zrv; %z reflections

    a(:,4:6) = P(i,4:6) + xyzrv;
    
    b = cat(1,b,a);
end
end


% function sa = fcnsolidangle(A,r)
% %A = area,  r = radius (same units as A)
% ct = 1; %cosine theta (1 for perpendicular)
% sa = 0.5*(1 - r./sqrt(r.^2+A/pi*ct));
% end
