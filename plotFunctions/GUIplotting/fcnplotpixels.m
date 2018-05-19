function fcnplotpixels(input,output,handles,G1,flags,PE,ip)

%PATCHES AND FACES --------------------------------------------------------
fcnpixels(input,output,handles,flags,PE,ip); %ip = photon indices

%PLOT SOLID ANGLES ------------------------------------------------------
delete(findobj(handles.GUI.figure1,'DisplayName','Projected Solid Angle'))
if flags.status.showSolidAngles
    fi = input.cube.all.fi; %face index
    for plane = 1:6
        if input.cube.pixelsperface(plane)>0
            colorMat = [    194 254 52
                110 41 247
                252 35 129
                0 0 204
                67 237 243
                215 213 17];
            color = colorMat(plane,:)/255;
            
            p0 = G1.p0;
            if input.cube.shapeID==3
                L = norm(-sign(p0+eps)*input.cube.height/2-p0);
            else
                L = 300;
            end
            
            n = 21;
            f = linspace(0,1,n)'; %fractions
            
            vx = input.cube.all.v.x(fi==plane,:)' - p0(1);  nx=size(vx,2);
            vy = input.cube.all.v.y(fi==plane,:)' - p0(2);
            vz = input.cube.all.v.z(fi==plane,:)' - p0(3);
            
            ov = ones(n,1); %onesvec
            vxu = [ov*vx(1,:)+f*(vx(2,:)-vx(1,:));  ov*vx(2,:)+f*(vx(3,:)-vx(2,:));  ov*vx(3,:)+f*(vx(4,:)-vx(3,:));  ov*vx(4,:)+f*(vx(1,:)-vx(4,:)); nan(1,nx)];
            vyu = [ov*vy(1,:)+f*(vy(2,:)-vy(1,:));  ov*vy(2,:)+f*(vy(3,:)-vy(2,:));  ov*vy(3,:)+f*(vy(4,:)-vy(3,:));  ov*vy(4,:)+f*(vy(1,:)-vy(4,:)); nan(1,nx)];
            vzu = [ov*vz(1,:)+f*(vz(2,:)-vz(1,:));  ov*vz(2,:)+f*(vz(3,:)-vz(2,:));  ov*vz(3,:)+f*(vz(4,:)-vz(3,:));  ov*vz(4,:)+f*(vz(1,:)-vz(4,:)); nan(1,nx)];
            [r,c] = size(vzu);
            
            sc = fcnCC2SCr(vxu(:),vyu(:),vzu(:));
            sc(:,1) = L;
            cc = fcnSC2CC(sc);
            x = reshape(cc(:,1),r,c)+p0(1);
            y = reshape(cc(:,2),r,c)+p0(2);
            z = reshape(cc(:,3),r,c)+p0(3);
            plot3(x(:),y(:),z(:),'-','color',color,'DisplayName','Projected Solid Angle');
            
%             %PLOT PIXEL 1 PYRAMID -----------------------------------------------------
%             p01 = ones(4,1)*G1.p0;
%             i = 1;
%             v1 = [1 11 22 33];
%             plot3([p01(:,1) x(v1,i)]',[p01(:,2) y(v1,i)]',[p01(:,3) z(v1,i)]','k');
            
            %     %PLOT PYRAMID SURFACE -----------------------------------------------------
            %     p01 = ones(5,1)*G1.p0;
            %     v1 = [1 11 22 33 1];
            %     h(5) = surf([p01(:,1) x(v1,i)],[p01(:,2) y(v1,i)],[p01(:,3) z(v1,i)],'FaceColor',[.7 .7 .7],'EdgeColor','none'); alpha(.8)
                        
            %      %PLOT SPHERE -----------------------------------------------------
            %      [x,y,z] = sphere(30); %
            %      x = x*L + p0(1);
            %      y = y*L + p0(2);
            %      z = z*L + p0(3);
            %      h(6) = surf(x,y,z,'FaceColor',[.6 .6 .6],'FaceAlpha',.1,'FaceLighting','gouraud','EdgeColor','None','AmbientStrength',1);
            %      %h = light('Position',[0 0 0],'Style','infinite'); handles.checkbox302 = [handles.checkbox302 h];
        end
    end
end

end


function fcnpixels(input,output,handles,flags,PE,ip)
if input.cube.pixels==0; return; end  
zm = zeros(input.cube.pixels,1);
Cdata = zm;  AlphaData = zm + .1;    

af = find(flags.status.plotside); %active faces
iap = find(any(input.cube.all.fi==af,2)); %active pixel indices

if flags.status.lappd
    if any(ip)
        nx = 100;
        sl = input.cube.pixelsize; %(mm) [length height]

        yv = linspace(-1,1,nx) * sl(1)/2;
        [c, strip, pid] = fcnanalogvoltage(input,flags,PE,input.plotTime,yv,ip,false);
        c = floorandceil(c,input.cube.dsp.dynamicrange);
        
        i = sum(c,2)>.05;  na = sum(i);
        if na>0
            c = c(i,:);
            pid = pid(i);            
            
            x = min(c,30)*0 + 1; %*0 flatten
            y = repmat(yv, [na 1]);
            z = repmat( input.cube.pmt.stripx(strip(i)), [1 nx]);

            rpy = repmat(fcnVEC2RPY(-input.cube.all.normalVec(pid,:)),[nx 1]);
            X = rotateW2Bc(rpy,[x(:) y(:) z(:)]);
            X = eshape(X,[na nx 3]) + reshape(input.cube.all.xyz(pid,:), [na 1 3]);   X(:,end+1,:)=nan;  c(:,end+1)=nan;
            x = X(:,:,1)';   y = X(:,:,2)';  z = X(:,:,3)';  c=c';
            
            k = c>.01; k([1 end-1 end],:)=true;  ov = [1 1];
            
            rpy = reshape(permute(reshape(rpy,[na nx 3]),[2 1 3]),[na*nx 3]);
            sw = fcnrotateW2B(rpy(:,1),rpy(:,2),rpy(:,3),[0 0 input.cube.pmt.stripwidth/2]);
            sw = reshape(sw,[nx numel(pid) 3]); sw(end+1,:,:)=nan;  sw=reshape(sw,(nx+1)*numel(pid),3);
            
            %h(1) = surface(x(k)*ov + sw(k,1)*[-1 1], y(k)*ov + sw(k,2)*[-1 1], z(k)*ov + sw(k,3)*[-1 1], c(k)*ov, 'facecolor','flat','edgecolor','none','Parent',handles.GUI.axes1);
            handles.stripanode.XData = x(k)*ov + sw(k,1)*[-1 1];
            handles.stripanode.YData = y(k)*ov + sw(k,2)*[-1 1];
            handles.stripanode.ZData = z(k)*ov + sw(k,3)*[-1 1];
            handles.stripanode.CData = c(k)*ov;
            AlphaData(pid) = 1;
        end
    end
    AlphaData = AlphaData*.25;
else
    if isempty(output(1).D); V=input.cube.dsp.N; else V=output(1).D; end
    if flags.status.plotcurrentvoltage
        Cdata = interp1(input.cube.dsp.t,V,input.plotTime)';
        AlphaData(abs(Cdata)>input.cube.dsp.noise*3) = .6;
    else
        %Cdata = sum(Vi,1)';
        Cdata = max(V,[],1)';
        %AlphaData(abs(Cdata)>min([inf; output(1).amplitude])) = .6;
         AlphaData(abs(Cdata)>.100) = .6; %>100 mV
    end
end
if any(iap)
    handles.gridanode.XData=input.cube.all.v.x(iap,:)';
    handles.gridanode.YData=input.cube.all.v.y(iap,:)';
    handles.gridanode.ZData=input.cube.all.v.z(iap,:)';
    handles.gridanode.CData=Cdata(iap,:);
    handles.gridanode.FaceVertexAlphaData=AlphaData(iap,:);
end

end