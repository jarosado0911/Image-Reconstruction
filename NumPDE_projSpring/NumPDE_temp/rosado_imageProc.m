function rosado_imageProc(noise,p,gamma,var)

if (nargin ==1)
    % Get Gaussion of initial raster image u0
    p=8.9; gamma=5.5*10^(-11.25);
    spd = 0.022; % adjusting speed of front
    var = 0.000015;
    %var = 0.000079;
end

%n = 200;       % number of space gridpoints
%[u_init,~,~]=u0(n,noise);

u_init = getImage();
[n,~]=size(u_init);

dt = 0.5e-3;    % time step
tf = 200e-0;    % final time
nr = 0;       % number of reinitialization steps
nsteps = 500;  % number of steps with graphic output
make_video = 1;
do_solve =1 ;
%-----------------------------------------------------------------------
% if nargin<1, method = 1; end
 nt = ceil(tf/dt); dt = tf/nt;
 x = linspace(0,1,n)'; h = x(2)-x(1);
 [X,Y] = meshgrid(x); ax = [min(x) max(x) min(x) max(x)];
%-----------------------------------------------------------------------

t = now;
myFolder = sprintf('Sim_Run_%s/Stills_%s', datestr(t,'mm-dd-yyyy HH-MM'),datestr(t,'mm-dd-yyyy HH-MM'));
if ~exist(myFolder, 'dir')
       mkdir(myFolder)
end
 
m=2;
% initial conditions
P = ((X-.5).^m+((Y-.5)).^m).^(1/m)*1.5-0.59;   

figure('units','normalized','outerposition',[0 0 0.5 0.5])
imagesc([0 1],[0 1]',u_init)
colorbar
xlabel('X');
ylabel('Y');
set(gca,'YDir','normal')
title(sprintf('Original Image, NxN = %i x %i',n,n));
filepath = sprintf('Sim_Run_%s/original.png', datestr(t,'mm-dd-yyyy HH-MM'));
saveas(gcf,filepath)

% % Get Gaussion of initial raster image u0
% p=8.0; gamma=5.5*10^(-11.25);
% spd = 0.052; % adjusting speed of front
% %var = 0.000065;
% var = 0.000025;
% %var = 0;

%Eu=gaussian_filter(u_init,X,Y',var);
Eu=gauss2d(u_init,x,n,var);

figure('units','normalized','outerposition',[0 0 0.5 0.5])
imagesc([0 1],[0 1],Eu)
colorbar
set(gca,'YDir','normal')
title(sprintf('Blurred, sigma = %0.2e, NxN = %i x %i',var,n,n))
xlabel('X');
ylabel('Y');
filepath = sprintf('Sim_Run_%s/blurred.png', datestr(t,'mm-dd-yyyy HH-MM'));
saveas(gcf,filepath)
%pause

[Eux, Euy]=gradient(Eu,h);

figure('units','normalized','outerposition',[0 0 0.75 0.75])
imagesc(x,x,sqrt(Eux.^2+Euy.^2));
colorbar
title(sprintf('|J*Du0|, var = %0.1e',var));
xlabel('X');
ylabel('Y');
set(gca,'YDir','normal')
filepath = sprintf('Sim_Run_%s/gradJDu0.png', datestr(t,'mm-dd-yyyy HH-MM'));
saveas(gcf,filepath)

% Osher/Fedkiw page 120
g = @(x,y) 1./(1+gamma*(abs(x).^2+abs(y).^2).^(p/2));  
figure('units','normalized','outerposition',[0 0 0.75 0.75])
imagesc(x,x,g(X,Y));
colorbar
caxis([0 1]);
set(gca,'YDir','normal')
title(sprintf('G(X,Y), p = %0.3f, gamma = %0.1e',p,gamma));
xlabel('X');
ylabel('Y');
filepath = sprintf('Sim_Run_%s/speedfunction_without_u0.png', datestr(t,'mm-dd-yyyy HH-MM'));
saveas(gcf,filepath)
%pause

g_DU = g(Eux,Euy);
A=g_DU;
figure('units','normalized','outerposition',[0 0 0.75 0.75])
imagesc([0 1], [0 1],A);
caxis([0 1]);
set(gca,'YDir','normal')
colorbar
title(sprintf('g(Du0), p = %0.3f, gamma = %f',p,gamma))
xlabel('X');
ylabel('Y');
filepath = sprintf('Sim_Run_%s/speedfunction_with_u0.png', datestr(t,'mm-dd-yyyy HH-MM'));
saveas(gcf,filepath)

%pause
%-----------------------------------------------------------------------
if (make_video == 1)
    fig=figure('units','normalized','outerposition',[0 0 1 0.95]);
    filepath = sprintf('Sim_Run_%s/rosado_levelset_vid', datestr(t,'mm-dd-yyyy HH-MM'));
    v = VideoWriter(filepath,'MPEG-4');
    open(v)
end

if do_solve == 1
count = 0;
for it = 1:nt 
   F = -curvature(A,P,h,spd)*5e-3;     % movement under curvature
   P = P-dt*FabsgradP(P,h,F);      % level set update
   
   % reinitialization steps
   if nr>0
    for ir = 1:nr                               
        P = P-dt*FabsgradP(P,h,P./sqrt(P.^2+(2*h)^2),1);
    end
   end
   
   % plotting purposes
   if it==1|floor(nsteps*it/nt)>floor(nsteps*(it-1)/nt) % visualization
        clf
        %subplot(1,2,1), 
        contourf(x,x,-P,[0 0],'k-')
        axis equal, axis(ax)
        title(sprintf('geometry at t=%0.2f',it*dt))
        
%         subplot(1,2,2), surf(x,x,-P,'EdgeAlpha',.2),
%         hold on
%         patch([0 1 1 0],[0 0 1 1],[0 0 0 0],'k','FaceAlpha',.5)
%         hold off, axp = [min(0,min(min(-P))) max(0,max(max(-P)))];
%         zlim([-1 0.6])
%         %axis([ax axp]), 
%         title('level set function')
        
        filepath = sprintf('Sim_Run_%s/Stills_%s/rosado_levelset_still_%i.png', datestr(t,'mm-dd-yyyy HH-MM'),datestr(t,'mm-dd-yyyy HH-MM'),count);
        saveas(gcf,filepath)
        count = count+1;
        if(make_video == 1)
            thisframe=getframe(fig);
            writeVideo(v, thisframe);
        end
        drawnow
   end
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

function F = curvature(A,P,h,spd)
% computes generalized curvature by central differences
Pxx = diff(P([1 1:end end],:),2)/h^2;
Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
Px = (P(3:end,:)-P(1:end-2,:))/(2*h); 
Px = Px([1 1:end end],:);
Py = (P(:,3:end)-P(:,1:end-2))/(2*h); 
Py = Py(:,[1 1:end end]);
Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); 
Pxy = Pxy(:,[1 1:end end]);

F1=(A.*((Pxx.*Py.^2-Px.*Py.*Pxy)+(Pyy.*Px.^2-Px.*Py.*Pxy)))./((Px.^2+Py.^2).^(1.5));

Ax = (A(3:end,:)-A(1:end-2,:))/(2*h); 
Ax = Ax([1 1:end end],:);
Ay = (A(:,3:end)-A(:,1:end-2))/(2*h); 
Ay = Ay(:,[1 1:end end]);

F2=(Ax.*Px+Ay.*Py)./((Px.^2+Py.^2).^(0.5));
F=(F1+F2);
F = min(max(F,-spd/h),spd/h);