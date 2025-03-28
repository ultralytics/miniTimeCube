% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function i = fcninsidedetector(x, input, ec)
%x is nx3 particle position
%ec is edge cut (optional)
if nargin==2;  ec=0;  end

if input.cube.shapeID==3 %cylinder
    i = sqrt(x(:,1).^2+x(:,2).^2)<=(input.cube.radius-ec) & abs(x(:,3))<=(input.cube.height/2-ec);
else %cube
    Lr = input.cube.Lr-ec;
    ax = abs(x);
    i = ax(:,1)<=Lr(1) & ax(:,2)<=Lr(2) & ax(:,3)<=Lr(3);
end
