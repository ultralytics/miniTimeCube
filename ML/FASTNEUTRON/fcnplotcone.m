% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [X,Y,Z,conevecs] = fcnplotcone(a,b,angle)
%function plots a cone with origin 'a', pointing toward 'b', half angle 'angle'
vec = b - a;
length = 0.6 * fcnrange(vec);
r = linspace(0,1*tand(angle),30);
C = fcnVEC2DCM_W2B(-vec);
C_M2NED = [ 0 1 0; 0 0 1; 1 0 0];%matlab cone to ned cone

[X,Y,Z] = cylinder(r,360);
sx = size(X);

c = [X(:) Y(:) Z(:)*1]*C_M2NED*C * length * sign(tand(angle));
X = reshape(c(:,1)+a(1),sx);
Y = reshape(c(:,2)+a(2),sx);
Z = reshape(c(:,3)+a(3),sx);

conevecs = zeros(sx(2),3);
conevecs(:,1) = X(end,:) - a(1);
conevecs(:,2) = Y(end,:) - a(2);
conevecs(:,3) = Z(end,:) - a(3);
%surf(X,Y,Z,'facecolor','r','edgecolor','none','facealpha',.3); axis equal vis3d; camlight headlight; camlight headlight
