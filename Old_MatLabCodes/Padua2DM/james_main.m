clear all; close all; more off;
clf
clc
myTestfun=@(x,y) (1+25.*x.^2.*y.^2).^(-1);

n=100;                 % PADUA POINTS DEGREE, for Lagrange interpolation.
xyrange=[-1,1,-1,1];   % DEFINITION OF THE RECTANGLE.

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

theM=[3:2:100];

for j=1:length(theM)
    clf
M=theM(j);                      
X = pdpts(M,xyrange);          
LnfX = pdint(n,xyrange,myTestfun,X);                         
fX=feval(myTestfun,X(:,1),X(:,2)); 


tri = delaunay(X(:,1),X(:,2));
light
subplot(1,2,1)
hold on
scatter(X(:,1),X(:,2),'.','r')
h=trisurf(tri, X(:,1),X(:,2),LnfX,'FaceLighting','gouraud','FaceColor','interp','AmbientStrength',0.5,'DiffuseStrength',1);
light('Position',[-1 -1 1],'Style','local')
set(h,'linestyle','none')
view(30,30);
subplot(1,2,2)
scatter(X(:,1),X(:,2),'.')

drawnow 
      %Capture the plot as an image 
      frame = getframe(hfig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
end