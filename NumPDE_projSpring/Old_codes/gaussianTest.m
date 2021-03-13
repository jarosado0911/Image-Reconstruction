function gaussianTest

npts=999;
[phi0,~,~,x,y]=initial_phi(npts,0,0);
dx = x(2)-x(1);
dy = y(2)-y(1);

[delphiX,delphiY,~,~]=gradPhi(npts,0);
[u_initial,~,~]=u0(npts);

[X,Y]=meshgrid(x,y);

E=gaussian_filter(phi0,X,Y',.015);
Ex=gaussian_filter(delphiX,X,Y',.015);
Ey=gaussian_filter(delphiY,X,Y',.015);
Eu=gaussian_filter(u_initial,X,Y',.015);

[DxJu, DyJu]=grad(Eu,dx,dy);

p=1;

normJu=zeros(npts,npts);

for j=1:npts
    for i=1:npts
        normJu(i,j)=sqrt(DxJu(j,i)^2+DyJu(j,i)^2)^p;
    end
end

gDelU =1./(1+normJu); 

[DelgDelUx, DelgDelUy]=grad(gDelU,dx,dy);


surf(X,Y,DelgDelUx)
view(0,90)
colormap jet
shading interp

% figure(1)
% surf(X,Y,Eu)
% view(0,90)
% colormap jet
% shading interp
% 
% figure(2)
% surf(X,Y,Ex)
% view(0,90)
% colormap jet
% shading interp
% 
% figure(3)
% surf(X,Y,Ey)
% view(0,90)
% colormap jet
% shading interp
% 
% figure(4)
% surf(X,Y,E)
% view(0,90)
% colormap jet
% shading interp

end

