function rosado_levelset_front(noise)
n = 120;       % number of space gridpoints
dt = 0.5e-3;    % time step
tf = 100e-0;    % final time
nr = 1;       % number of reinitialization steps
nsteps = 500;  % number of steps with graphic output
%-----------------------------------------------------------------------
% if nargin<1, method = 1; end
 nt = ceil(tf/dt); dt = tf/nt;
 x = linspace(0,1,n)'; h = x(2)-x(1);
 [X,Y] = meshgrid(x); ax = [min(x) max(x) min(x) max(x)];
%-----------------------------------------------------------------------
% initial conditions
P = sqrt((X-.5).^2+((Y-.5)).^2)*1.5-0.45;   

[u_init,~,~]=u0(n,noise);
figure('units','normalized','outerposition',[0 0 0.75 0.55])
subplot(1,2,1)
imagesc([0 1],[0 1]',u_init)
set(gca,'YDir','normal')
title('Original Image')
% Get Gaussion of initial raster image u0
%var = 0.0075*10^(0.6); %0.015*10^(0.6);
var = 0.05;
Eu=gaussian_filter(u_init,X,Y',var);
subplot(1,2,2)
imagesc([0 1],[0 1],Eu)
set(gca,'YDir','normal')
title(sprintf('Blurred, sigma = %0.3f',var))
pause

% Osher/Fedkiw page 120
[Eux, Euy]=gradient(Eu);
% Eux = (Eu(3:end,:)-Eu(1:end-2,:))/(2*h); 
% Eux = Eux([1 1:end end],:);
% Euy = (Eu(:,3:end)-Eu(:,1:end-2))/(2*h); 
% Euy = Euy(:,[1 1:end end]);
Eux = min(max(Eux,-1/h),1/h);
Euy = min(max(Euy,-1/h),1/h);

%EuxNorm = Eux;%./max(max(Eux));
%EuyNorm = Euy;%./max(max(Euy));

% figure('units','normalized','outerposition',[0 0 1 0.5])
% subplot(1,3,1)
% imagesc([0 1],[0 1], abs(EuxNorm));
% set(gca,'YDir','normal')
% title('D(u_{blurred})_x')
% colorbar
% subplot(1,3,2)
% imagesc([0 1],[0 1], abs(EuyNorm));
% set(gca,'YDir','normal')
% title('D(u_{blurred})_y')
% colorbar
% subplot(1,3,3)
% imagesc([0 1],[0 1], (EuxNorm.^2+EuyNorm.^2).^0.5);
% set(gca,'YDir','normal')
% title('||Du_{blurred}||')
% colorbar
% pause

%Enorm =(EuxNorm.^2+EuyNorm.^2).^0.5;

%Eux(Enorm<3.5e-1)=0;
%Euy(Enorm<3.5e-1)=0;

% figure
% subplot(1,3,1)
% imagesc([0 1], [0 1], abs(Eux));
% colorbar
% subplot(1,3,2)
% imagesc([0 1], [0 1], abs(Euy));
% colorbar
% subplot(1,3,3)
% imagesc([0 1], [0 1], (Eux.^2+Euy.^2).^0.5);
% colorbar

%pause
p=2.5001;%0.5;
alpha=1;
beta=1;
gamma=5;
p1 = 1000000000;
p2 = 1000000000;

g = @(x,y) alpha./(beta+gamma*(abs(x).^p1+abs(y).^p2).^p);  
figure
imagesc(x,x,g(X,Y));
pause
g_DU = g(Eux,Euy);
%g_DU=1./(0.5+3*(Eux.^4+Euy.^4).^p);

A=g_DU;
figure('units','normalized','outerposition',[0 0 0.4 0.5])
imagesc([0 1], [0 1],A);
colorbar
title('g(Du0)')
set(gca,'YDir','normal')
pause

%A(A>0.65)=3;
% figure('units','normalized','outerposition',[0 0 0.4 0.5])
% imagesc([0 1], [0 1],A);
% colorbar
% title('g(Du0)')
% set(gca,'YDir','normal')
% pause

% [Ax,Ay]=gradient(A);
% figure('units','normalized','outerposition',[0 0 1 0.5])
% subplot(1,3,1)
% imagesc([0 1], [0 1], abs(Ax));
% title('Speed in x')
% set(gca,'YDir','normal')
% colorbar
% subplot(1,3,2)
% imagesc([0 1], [0 1],abs( Ay));
% title('Speed in y')
% set(gca,'YDir','normal')
% colorbar
% subplot(1,3,3)
% imagesc([0 1], [0 1], (Ax.^2+Ay.^2).^0.5);
% title('||Speed||')
% set(gca,'YDir','normal')
% colorbar
% 
% pause
%-----------------------------------------------------------------------
fig=figure('units','normalized','outerposition',[0 0 1 0.95]);
v = VideoWriter('rosado_levelset_test5','MPEG-4');
open(v)
for it = 1:nt 
   F = -curvature(A,P,h)*5e-3;     % movement under curvature
   P = P-dt*FabsgradP(P,h,F);                        % level set update
%     for ir = 1:nr                               % reinitialization steps
%         P = P-dt*FabsgradP(P,h,P./sqrt(P.^2+(2*h)^2),1);
%     end

    [Px, Py] = gradient(P);
    
   if it==1|floor(nsteps*it/nt)>floor(nsteps*(it-1)/nt) % visualization
        clf
        subplot(1,2,1), 
        contourf(x,x,-P,[0 0],'k-')
        axis equal, axis(ax)
        title(sprintf('geometry at t=%0.2f',it*dt))
        
        subplot(1,2,2), surf(x,x,-P,'EdgeAlpha',.2),
        hold on
        patch([0 1 1 0],[0 0 1 1],[0 0 0 0],'k','FaceAlpha',.5)
        hold off, axp = [min(0,min(min(-P))) max(0,max(max(-P)))];
        zlim([-0.3 0.3])
        axis([ax axp]), title('level set function')
        
%         subplot(1,3,3)
%         quiver(x,x, (Px)',(Py)',0)
        
        thisframe=getframe(fig);
        writeVideo(v, thisframe);
        drawnow
   end
end

%=======================================================================

function dP = FabsgradP(P,h,F,c)
if nargin<4, c = 0; if nargin<3, F = 1; end, end
DxP = diff(P)/h;   DxmP = DxP([1 1:end],:); DxpP = DxP([1:end end],:);
DyP = diff(P')'/h; DymP = DyP(:,[1 1:end]); DypP = DyP(:,[1:end end]);
Np = sqrt(max(DxmP,0).^2+min(DxpP,0).^2+max(DymP,0).^2+min(DypP,0).^2);
Nm = sqrt(min(DxmP,0).^2+max(DxpP,0).^2+min(DymP,0).^2+max(DypP,0).^2);
dP = max(F,0).*(Np-c)+min(F,0).*(Nm-c);

function F = curvature(A,P,h)
% computes curvature by central differences
Pxx = diff(P([1 1:end end],:),2)/h^2;
Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); 
Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); 
Py = Py(:,[1 1:end end]);
Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); 
Pxy = Pxy(:,[1 1:end end]);

%F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;

F1=(A.*((Pxx.*Py.^2-Px.*Py.*Pxy)+(Pyy.*Px.^2-Px.*Py.*Pxy)))./((Px.^2+Py.^2).^(1.5));

Ax = (A(3:end,:)-A(1:end-2,:))/(2*h); 
Ax = Ax([1 1:end end],:);
Ay = (A(:,3:end)-A(:,1:end-2))/(2*h); 
Ay = Ay(:,[1 1:end end]);

F2=(Ax.*Px+Ay.*Py)./((Px.^2+Py.^2).^(0.5));
beta = 1; % = 0.5
F=(F1+F2)*beta;
%F = min(max(F,-1),1);
F = min(max(F,-0.5/(h)),0.5/(h));