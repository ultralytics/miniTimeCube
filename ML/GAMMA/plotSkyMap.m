function [] = plotSkyMap(F,v,angles,np,particleName)
nr = 100;
nc = 200;
r = linspace(90,-90,nr);
c = linspace(-180,180,nc); 
[rm,cm] = meshgrid(r,c);
cc = fcnSC2CCd(1,rm(:),cm(:));
z = zeros(nr*nc,1);

particleName=strrep(particleName,'Fiber','');
coneFlag=any(angles>0);

n=size(v,1);
v=fcnvec2uvec(v);
for i=1:min(n,np)
    ct=sum(cc.*v(i,:),2);
    if coneFlag
        ct=cos( abs(acos(ct)-d2r*angles(i)) );
    end
    z = z + F(ct).^4;
end
zm = reshape(z,nc,nr);
%SUBTRACT APRIORI BACKGROUND MAP
%save zmBACKGROUND.mat zm
B=load('zmBACKGROUND.mat');  %size 200x100
zm = zm/sum3(zm);
B.zm=B.zm/sum3(B.zm);
zm=zm-B.zm;


%PLOT ---------------------------------------------------------------------
sca; sca; 
%fcnPlateCaree2sphere(rm,cm,zm)

%NEW CODE
sca;
[rmp,cmp] = fcnmapprojection(rm',cm','winkeltripel');
pcolor(cmp,rmp,zm'); shading flat;
input=[]; input = evalin('caller','input');
fcnplotdetectorprojection(input,[],ones(1536,1),[0 0 0],'winkeltripel','PMT');  %title(sprintf('%g %s',ne,s));
h=gca; h.CameraViewAngle=7.5; h.Position(1)=h.Position(1)-.02;
hp=findobj(gca,'Type','Patch'); hp.FaceAlpha=.3;
fcnfontsize(12); h.CLim=minmax3(zm);
%NEW CODE


%OLD CODE
%pcolor(cm,rm,zm); 
%plot(0,0,'+','MarkerSize',50,'color',[.7 .7 .7]*0,'linewidth',2);
%shading flat; axis equal tight off
%s=''; if min(n,np)>1; s='s'; end;  text(-179,89,sprintf('%g %s%s',min(n,np),particleName,s),'Color',[1 1 1],'Fontsize',54,'VerticalAlignment','Top')
%OLD CODE



switch particleName
    case 'Gamma'
        colormap(hot);
    case 'Neutron'
        colormap(parula);
end

end

function fcnPlateCaree2sphere(rm,cm,zm)
sr = 1;
np1 = 600;
np2 = np1+1;
[sx,sy,sz] = sphere(np1);
sc = fcnCC2SC(sx(:),sy(:),sz(:));
el2 = sc(:,2); az2 = sc(:,3);
scL = interp2(rm,cm,zm,reshape(el2,np2,np2),reshape(az2,np2,np2),'linear');
scL = reshape(scL,np2,np2);
sx=sx*sr;
sy=sy*sr;
sz=sz*sr;
surf(sx,sy,sz,scL); shading flat; caxis([min3(zm) max3(zm)+eps])
set(gca,'YDir','Reverse','ZDir','Reverse');
view(37,30)

fcnplotaxes3(sr*1.25)
fcnplotspherecross(0,0,30,[0 0 0],sr*1.02)
axis off tight equal vis3d
end