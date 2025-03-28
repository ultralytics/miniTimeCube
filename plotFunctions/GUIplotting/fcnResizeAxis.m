% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnResizeAxis(flags,handles,input,G1)
h = handles.GUI.axes1;
padding = [-1 1 -1 1 -1 1]*3; %(mm) extra empty border around event

if get(handles.GUI.checkboxautoresize,'Value')==1
    %a = axis; %current axis
    uvx = get(handles.uvec(1),{'XData','YData','ZData'}); set(handles.uvec(1),{'XData','YData','ZData'},{0,0,0});
    a = cell2mat(get(h,{'xlim','ylim','zlim'}));
    [mmx,mmy,mmz] = fcnaxesdatalims(h);    
    b = [mmx, mmy mmz]; %desired axis
    set(handles.uvec(1),{'XData','YData','ZData'},uvx);
    
    b(1) = min(b(1),min(uvx{1}));
    b(2) = max(b(2),max(uvx{1})); 
    b(3) = min(b(3),min(uvx{2})); 
    b(4) = max(b(4),max(uvx{2})); 
    b(5) = min(b(5),min(uvx{3})); 
    b(6) = max(b(6),max(uvx{3})); 

    if any(a~=(b+padding))
        %RESIZE UVEC0    
        c = [mmx; mmy; mmz];
        p0 = G1.p0;
        
        i=[1 1; 2 2; 3 3];
        uvec = -G1.uvec(1,i);
        d = ((c+ones(3,1)*[-1 1]) - p0(i))./uvec(i);   d(d<0)=inf;
        [minval, j] = min3(d);  
        
        extension = diff(c(j(1),:))*.15; %.15 extension past event
        p = [p0; p0-G1.uvec(1,:)*max(minval+extension,10)]; %10mm min length
        set(handles.uvec(1),'XData',p(:,1),'YData',p(:,2),'ZData',p(:,3))

        
        %FIND OUT HOW FAST TO ZOOM IN OR OUT
        [mmx,mmy,mmz] = fcnaxesdatalims(h);    b = [mmx mmy mmz];
        d0 = abs([a(1)-a(2) a(3)-a(4) a(5)-a(6)]);
        d = abs(a-b); d = [d(1)+d(2) d(3)+d(4) d(5)+d(6)];
        maxPercent=max(d./d0)*100; %maximum percent change on any axis
        np = max(min(maxPercent, 12),6);
        np = 1;
    
        %ZOOM IN OR OUT
        %a = a + padding; 
        b = b + padding;
        for n = [logspace(0,-2,np) 0] %interpolate between the two "zoom in" or "zoom out"
            axis(h, a*(n)+b*(1-n)); 
            pause(.01)
        end
    end

end
