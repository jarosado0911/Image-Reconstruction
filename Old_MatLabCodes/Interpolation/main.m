%REGULAR INTERPOLATION
clf
close all
clc
clear
myfun=@(x) (1+25*x.^2).^(-1);
a=-1;
b=1;
nsamples=(1:2:40);
frm=ceil(length(nsamples)/2);

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

choose=input('Equidistant[1] or Tchebyshev[2]: ');
for k=1:length(nsamples)
    clf
Np1=nsamples(k)
h=abs(b-a)/(Np1-1);
xlst=linspace(a,b,Np1);
flst=myfun(xlst);

clst=newtonCoeff(xlst,flst);
pxlst=linspace(a,b,1000);
pylst=newtonP(pxlst,xlst,clst);

xTcheb=(1/2)*(a+b)+(1/2)*(b-a)*cos( (2* (0:Np1-1)+1) / (2*Np1) * pi );
fTcheb=myfun(xTcheb);
cTcheb=newtonCoeff(xTcheb,fTcheb);
yTcheb=newtonP(pxlst,xTcheb,cTcheb);

if choose == 2
%subplot(2,frm,k)
hold on
scatter(xTcheb,fTcheb,'filled','b')
plot(pxlst,yTcheb,'r')
plot(pxlst,myfun(pxlst),'Black')
ylim([0 1.25])

drawnow 
     % Capture the plot as an image 
      frame = getframe(hfig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if k == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
else
    
%subplot(2,frm,k)
hold on 
scatter(xlst,flst,'filled','b')
plot(pxlst,pylst,'r') 
plot(pxlst,myfun(pxlst),'Black')
ylim([0 1.25])

drawnow 
     % Capture the plot as an image 
      frame = getframe(hfig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if k == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
end

%QUESTION #3
%For the interpolation using equidistant points, it looks like that
%convergence will not occur near and at the endpoints.

%For the interpolation using Tchebyshev points, convergence is much better,
%but I noticed that for more points, there is some disparity at x=-1.