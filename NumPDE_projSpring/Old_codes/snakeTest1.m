function snakeTest1()

npts = 599;

% Initial image
[u_initial,x,y]=u0(npts);

% mesh grid for surf plots and Gaussian Blurring
[X,Y]=meshgrid(x,y);

    surf(X,Y,u_initial)
    view(0,90)
    colormap jet
    shading interp
    pause


% grid resolutions
dx = x(2)-x(1);
dy = y(2)-y(1);

% Get Gaussion of initial raster image u0
var = 0.03; %0.28; %0.5;
Eu=gaussian_filter(u_initial,X,Y',var);

surf(X,Y,Eu)
    view(0,90)
    colormap jet
    shading interp
    pause

% This is for the edge detection coefficients as described 
% Osher/Fedkiw page 120
[DxJu, DyJu]=grad(Eu,dx,dy);
p=1;

normJu=zeros(npts,npts);
for j=1:npts
    for i=1:npts
        normJu(i,j)=sqrt(DxJu(j,i)^2+DyJu(j,i)^2)^p;
    end
end

g_DU =1./(1+normJu); 

% This is for equation 12.3 on page 120
[Dxg_DU,Dyg_DU]=grad(g_DU,dx,dy);

Dxg_DU=Dxg_DU/200;
Dyg_DU=Dyg_DU/200;

%pause
%Dxg_DU=20*abs(Dxg_DU);
%Dyg_DU=20*abs(Dyg_DU);

% surf(X,Y,Dxg_DU)
% view(0,90)
% colormap jet
% shading interp
% pause

% initial curve
r0 = 0.5;
[phi0,~,~,~,~]=initial_phi(npts,r0,0);
%[phi0,~,~,~,~]=initial_sqr(npts,r0,0);

k=0.9*dx;
cfl=k/dx;

tf = 5;
nT = floor(tf/k);

% initial curve
phi = phi0;
phi_old=phi;

figure
% advection solve only

fig=figure;
v = VideoWriter('my_first_snake2','MPEG-4');
open(v)

for n=1:nT

    for i=2:npts-1
    
        for j=2:npts-1
        
            phi(i,j) = phi(i,j)+cfl*(Dxg_DU(i,j)*(phi(i,j+1)-phi(i,j-1))...
                                    +Dyg_DU(i,j)*(phi(i+1,j)-phi(i-1,j)));
                         
            if phi(i,j)>10
                phi(i,j)=1;
            end
        end
    end
   max(max(abs(phi-phi_old)))
    if max(max(abs(phi-phi_old)))<=1.06%1.4
        phi_old = phi;
        r0 = r0*0.99;
        [phi,~,~,~,~] = initial_phi(npts,r0,0);
    
    end
    
    surf(X,Y,phi)
    title(sprintf('t = %0.4f',n*k));
    view(0,90)
    colormap jet
    shading interp
    
    thisframe=getframe(fig);
    writeVideo(v, thisframe);
    drawnow
end

end

