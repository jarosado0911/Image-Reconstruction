function [delphiX,delphiY,x,y]=gradPhi(npts,grph)

%npts = 1200;

[phi0,~,~,x,y]=initial_phi(npts,grph);

nX = length(x);
nY = length(y);

hx = x(2)-x(1);
hy = y(2)-y(1);

delphiX = zeros(nX,nY);
delphiY = zeros(nX,nY);

normDelphi=zeros(nX,nY);

for i=2:nX-1
    for j=2:nY-1
        delphiX(i,j) = (phi0(i+1,j)-phi0(i-1,j))/(2*hx);
        delphiY(i,j) = (phi0(i,j+1)-phi0(i,j-1))/(2*hy);
        normDelphi(i,j)=norm([delphiX(i,j) delphiY(i,j)],2);
    end
end

% figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.75]);
% [X,Y]=meshgrid(x(2:nX-1),y(2:nY-1));
% surf(X,Y,normDelphi(2:nX-1,2:nY-1));
% xlim([-0.025 1.025]);
% ylim([-0.025 1.025]);
% view(0,90)
% colormap jet
% shading interp
end

