function levelSetv2()

npts = 299;
% Initial image
[u_initial,x,y]=u0(npts);
% RGB = imread('grumpycat.png');
% 
% u_initial = double(rgb2gray(RGB)./200);
% [height, width] = size(u_initial);
% 
% x=[1:width];
% y=[1:height];

% mesh grid for surf plots and Gaussian Blurrin
[X, Y] = meshgrid(x, y);
nptsX = length(x); % height;
nptsY = length(y); % width;

% grid resolutions
dx = 1/(x(2)-x(1));
%dy = y(2)-y(1);

% Get Gaussion of initial raster image u0
var = 0.1e2;
Eu=gaussian_filter(u_initial,X,Y',var);

% This is for the edge detection coefficients as described 
% Osher/Fedkiw page 120
[DxJu, DyJu]=gradient(Eu);
p=5;

%for coloring purposes
c = [255, 130, 255]*1/255;

normJu=zeros(nptsX,nptsY);
for i=1:nptsX
    for j=1:nptsY
        normJu(i,j)=sqrt(DxJu(i,j)^2+DyJu(i,j)^2)^p;
    end
end

g_DU =1./(1+normJu); 

% This is for equation 12.3 on page 120
[Dxg_DU,Dyg_DU]=gradient(g_DU);

% initial level set function
%phi0 = sqrt((X-width/2).^2+(Y-height/2).^2)-120*2;
phi0 = sqrt((X-0.5).^2+(Y-0.5).^2)-1;
%phi0 = 9*X+2*Y-0.5;

k=0.00095;
cfl=k/dx;
cfl
k
tf = 2.0;
nT = floor(tf/k);

% initial curve
phi = phi0;
phi2 = phi0;

%fig=figure('units','normalized','outerposition',[0 0 0.65 0.475]);
fig=figure('units','normalized','outerposition',[0 0 0.75 0.95]);
v = VideoWriter('my_first_levelSet_testb','MPEG-4');
open(v)
kurv=k*curvature(u_initial,dx);

for n=1:nT
    %[d1,d2]=gradient(phi);
    for i=2:nptsX-1
        for j=2:nptsY-1
            phi(i,j)=cfl*sqrt((phi(i,j)-phi(i,j-1))^2+(phi(i,j)-phi(i-1,j))^2)*g_DU(i,j)*kurv(i,j)...
                +(1+cfl*(Dxg_DU(i,j)+Dyg_DU(i,j)))*phi(i,j)-cfl*(Dxg_DU(i,j)*phi(i-1,j)+Dyg_DU(i,j)*phi(i,j-1));
            
            phi2(i,j)=(1+cfl*(Dxg_DU(i,j)+Dyg_DU(i,j)))*phi2(i,j)-cfl*(Dxg_DU(i,j)*phi2(i-1,j)+Dyg_DU(i,j)*phi2(i,j-1));
        end
    end
        subplot(2,2,1)
        imagesc(x,y,u_initial)
        title('Original Image')
        colormap jet
        shading interp
    
        subplot(2,2,2)
        imagesc(x,y,Eu)
        title(sprintf('Gaussian Blurred var = %0.4f',var))
        colormap jet
        shading interp

        subplot(2,2,3)
        %imagesc(x,flip(y),phi)
        %axis xy
        %hold on
        contour(x,flip(y),phi,[-10 -10],'k-')
        %hold off    
        caxis([-100,100]*.5)
        title(sprintf('Time evolution with Curvature t = %0.4f',n*k));
        colormap jet
        shading interp

        subplot(2,2,4)
        imagesc(x,flip(y),phi2)
        axis xy
        hold on
        contour(x,flip(y),phi2,[0 0],'color',c,'LineWidth',0.95)
        hold off    
        caxis([-100,100]*.5)
        title(sprintf('Time evolution noCurvature t = %0.4f',n*k));
        colormap jet
        shading interp
        
%         if n==1
%             saveas(fig,'figt1','png');
%         end
%         
%         if n==33
%             saveas(fig,'figt33','png');
%         end
%         
%         if n==66
%             saveas(fig,'figt66','png');
%         end
%         
%         if n==102
%             saveas(fig,'figt102','png');
%         end
        
        thisframe=getframe(fig);
        writeVideo(v, thisframe);
        
        drawnow
end



end

function F = curvature(P,h)
% computes curvature by central differences
Pxx = diff(P([1 1:end end],:),2)/h^2;
Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); Py = Py(:,[1 1:end end]);
Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); Pxy = Pxy(:,[1 1:end end]);
F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;
F = min(max(F,-1/h),1/h);
end
