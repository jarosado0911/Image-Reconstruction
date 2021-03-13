function rosado_levelset_front(method)
%MIT18086_LEVELSET_FRONT
%    Computes the movement of fronts under a given velocity field
%    using the level set method. The iteration scheme is first order
%    accurate, and a reinitialization iteration is used.
%    The front velocity can be chosen to be constant or proportional
%    to the negative curvature.

% 03/2008 by Benjamin Seibold
% Feel free to modify for teaching and learning.
%-----------------------------------------------------------------------
n = 300;       % number of space gridpoints
dt = 1e-3;    % time step
tf = 100e-0;    % final time
nr = 1;       % number of reinitialization steps
nsteps = 100;  % number of steps with graphic output
%-----------------------------------------------------------------------
if nargin<1, method = 1; end
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,1,n)'; h = x(2)-x(1);
[X,Y] = meshgrid(x); ax = [min(x) max(x) min(x) max(x)];
%-----------------------------------------------------------------------
% initial conditions
P = sqrt((X-.5).^2+((Y-.5)).^2)-0.5;         
%P = 0.5*sin(2*pi*X).*sin(2*pi*Y)-0.25;
%P = max(P,.07-sqrt((X-.6).^2+(Y-.78).^2));     % cut a hole out
%P = min(P,max(max(.7-X,X-.87),max(.2-Y,Y-.6))); % add a rectangle
%P = min(P,sqrt((X-.3).^2+(Y-.25).^2)-.195);    % add a circle

A = P*0+1; % Make field of ones
[u_init,~,~]=u0(n);
figure
imagesc([0 1],[0 1],u_init)
% Get Gaussion of initial raster image u0
var = 0.001e1;
Eu=gaussian_filter(u_init,X,Y',var);
figure
imagesc([0 1],[0 1],Eu)
pause

% Osher/Fedkiw page 120
[DxJu, DyJu]=gradient(Eu);
p=5;

normJu=zeros(n,n);
for i=1:n
    for j=1:n
        normJu(i,j)=sqrt(DxJu(i,j)^2+DyJu(i,j)^2)^p;
    end
end

g_DU =1./(1+normJu); 

A=g_DU;
%-----------------------------------------------------------------------
fig=figure('units','normalized','outerposition',[0 0 0.75 0.95]);
v = VideoWriter('rosado_levelsetfront','MPEG-4');
open(v)
for it = 1:nt
    switch method
        case 1, F = 5e-2;              % movement with constant velocity
        case 2, F = -curvature(A*0+1,A+P,h)*5e-3;     % movement under curvature
    end
    P = P-dt*FabsgradP(P,h,F);                        % level set update
    for ir = 1:nr                               % reinitialization steps
        P = P-dt*FabsgradP(P,h,P./sqrt(P.^2+(2*h)^2),1);
    end
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
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); Py = Py(:,[1 1:end end]);
Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); Pxy = Pxy(:,[1 1:end end]);
%F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;

F=(A.*(Pxx.*Py.^2-Px.*Py.*Pxy)+(Pyy.*Px.^2-Px.*Py.*Pxy))./(Px.^2+Py.^2).^1.5;

F = min(max(F,-1/h),1/h);