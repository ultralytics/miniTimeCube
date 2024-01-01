function fcnplotdetectorprojection(input,iap,fxva,p0,proj,label)
if ~exist('iap','var')   || isempty(iap);     iap = 1:input.cube.pixels;      end %iap = indices of active pixels
if ~exist('fxva','var')  || isempty(fxva);    fxva = ones(size(iap))*100;     end %fxva = values of pixels
if ~exist('p0','var')    || isempty(p0);      p0 = [0 0 0];                   end %p0 = vertex for projection
if ~exist('proj','var')  || isempty(proj);    proj = 'winkeltripel';          end
if ~exist('label','var') || isempty(label);   label = 0;                      end %label S and PMT flag
fxva=fxva(:); %iap=iap(:);

interpflag = false; %interpolate dead pixels from neighbors
if interpflag
    j=isnan(fxva);  i=~j; %#ok<UNRCH>
    if any(j) && sum(i)>15 %at least 15 good pixels
        a=input.cube.all.xyz;
        F=scatteredInterpolant(a(i,1),a(i,2),a(i,3),fxva(i),'natural','linear');
        fxva(j) = F(a(j,1),a(j,2),a(j,3));
    end
    %iap=true(size(fxva));
end


x = input.cube.all.v.x - p0(1);
y = input.cube.all.v.y - p0(2);
z = input.cube.all.v.z - p0(3);
n = numel(fxva);

sc = fcnCC2SC(x(:),y(:),z(:));
rm = reshape(sc(:,2),[n 4]);
cm = reshape(sc(:,3),[n 4]); 

i=find(cm(:,1)<-100 & cm(:,2)>100);  %j = cm(i,1)==-180;  j=i(j);  i=i(~j);  cm(j,[1 4]) = cm(j,[1 4])+360;
if any(i) %pixels wrapping around the 180deg line
    if all(p0==0) %we can split pixels, but p0 must be in center of detector to make left and right symmetric!
        cm(i,2:3)=-180;  cm=[cm; -cm(i,:)];  rm=[rm; rm(i,:)];  fxva=[fxva; fxva(i)]; %SPLIT PIXELS IN HALF LEFT AND RIGHT
    else
        cm(i,2:3) = cm(i,2:3)-360; %EXTEND WRAPPED PIXELS TO THE LEFT
    end
end

[rm,cm] = fcnmapprojection(rm,cm,proj);  i=isnan(fxva); ad=ones(size(fxva)); ad(i)=.10; fxva(i)=min(fxva);
%if ~isempty(input.MTC);  i=fcnpruninglist; rm=rm(i,:); cm=cm(i,:); fxva=fxva(i); ad=ad(i); end %REMOVE PRUNED PIXELS FROM PLOT
if ~isempty(input.MTC);  i=fcnpruninglist; ad(:)=.1; ad(i)=1; end %DIM PRUNED PIXELS
patch(cm', rm', rm'*0,'w', ...
        'EdgeColor','none', ...
        'EdgeAlpha',0, ...
        'FaceColor','Flat', ...
        'FaceAlpha','Flat', ... %.2
        'FaceVertexCData',fxva, ...
        'FaceVertexAlphaData',ad, ...
        'CDataMapping','scaled', ...
        'AlphaDataMapping','none', ...
        'FaceLighting','none', ...
        'BackFaceLighting','unlit');

% for i=1:1536
%    text(mean(cm(i,:)),mean(rm(i,:)),sprintf('%g',i),'color',[1 1 1]*.3,'fontsize',12,'horizontalalignment','center','verticalalignment','middle');
% end
    
%PLOT BORDER
np = 101;
ov = ones(1,np); %onesvec
vx = 180*[-1*ov   linspace(-1,1,np)   1*ov   linspace(1,-1,np)]';
vy =  90*[linspace(-1,1,np)   1*ov   linspace(1,-1,np)  -1*ov ]';
[rm,cm] = fcnmapprojection(vy,vx,proj);
plot(cm,rm,'-','color',[1 1 1]*.7,'MarkerSize',2);

box = input.cube.box;  n=numel(box.x);  C=cell(3,n);
box.x = box.x - p0(1);
box.y = box.y - p0(2);
box.z = box.z - p0(3);
for i=1:n-1
    if isnan(box.x(i));
        C(:,i)={nan};
    else
        %r = rangec([box.x(i) box.y(i) box.z(i)],[box.x(i+1) box.y(i+1) box.z(i+1)]);
        C{1,i} = linspace(box.x(i),box.x(i+1),60);
        C{2,i} = linspace(box.y(i),box.y(i+1),60);
        C{3,i} = linspace(box.z(i),box.z(i+1),60);
    end
end
[el,az]=fcnelaz(cat(2,C{1,:}), cat(2,C{2,:}), cat(2,C{3,:}));
i=[false abs(diff(az)) > pi];  az(i)=nan; %wrapping values
[rm,cm] = fcnmapprojection(el*r2d,az*r2d,proj);
plot(cm(:),rm(:),'-','color',[1 1 1]*.7,'linewidth',1);

h=gca;
fcntight(h,'xy csigma');
daspect([1 1 1])
h.Visible='off'; 
h.Title.Visible='on';
%axis off tight equal

%LABEL SCROD AND PMT
if any(label)
    try
        switch label
            case 'SCROD'; vi=1;
            case 'ASIC row'; vi=2; return
            case 'ASIC column'; vi=3; return
            case 'ASIC channel'; vi=4; return
            case {'PMT','MCP'}; vi=5;
            case 'ASIC'; vi=6; return
            case 'Detector ASIC'; vi=7;
            case 'channel'; vi=8;
        end
        %vi = [1 5];
        labels = {'S','Ar','Ac','Ch','P','A','DA',''};
        labels = labels(vi); %#ok<*NASGU>
        
        ni = numel(vi);
        SRCCHi = MTCpixelID2SRCCH(1:1536);
        for i=1:ni %SCROD AND PMT
            si = SRCCHi(:,vi(i));  if vi==8; si=1:input.cube.pixels; end
            ui = fcnrow(fcnunique(si));  nui=numel(ui);  x=zeros(nui,3); s=cell(nui,1);
            for j=1:nui
                x(j,:) = mean(input.cube.all.xyz(si==ui(j),:),1);
                %s{j} = sprintf('%s%g',labels{i},ui(j));
                s{j} = sprintf('%s%g','',ui(j));
            end
            sc = fcnCC2SC(x(:,1)-p0(1),x(:,2)-.001-p0(2),x(:,3)-p0(3));
            [rm,cm] = fcnmapprojection(sc(:,2),sc(:,3),proj);
            text(cm,rm,s,'color',[1 1 1]*.3,'fontsize',18,'horizontalalignment','center','verticalalignment','middle');
        end
    end
end

end
