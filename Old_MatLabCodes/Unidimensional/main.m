%REGULAR INTERPOLATION
clf
close all
clc
clear
my2dfun=@(x,y) (1+25.*x.^2.*y.^2).^(-1);
a=-1;
b=1;
%nsamples=(1:2:40);
nsamples=24;
frm=ceil(length(nsamples)/2);

yslice=linspace(a,b,nsamples);

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

choose=input('Equidistant[1] or Tchebyshev[2]: ');
for k=1:length(yslice)
    %clf
    myfun=@(x) my2dfun(x,yslice(k));
%Np1=nsamples(k);
Np1=nsamples;
h=abs(b-a)/(Np1-1);
xlst=linspace(a,b,Np1);
flst=myfun(xlst);
yslice(k)
clst=newtonCoeff(xlst,flst);
pxlst=linspace(a,b,1000);
pylst=newtonP(pxlst,xlst,clst);

xTcheb=(1/2)*(a+b)+(1/2)*(b-a)*cos( (2* (0:Np1-1)+1) / (2*Np1) * pi );
fTcheb=myfun(xTcheb);
cTcheb=newtonCoeff(xTcheb,fTcheb);
yTcheb=newtonP(pxlst,xTcheb,cTcheb);

yslicelst=yslice(k)*ones(1,length(xlst));
pyslicelst=yslice(k)*ones(1,length(pxlst));

if choose == 2
%subplot(2,frm,k)
hold on
scatter3(yslicelst,xTcheb,fTcheb,'filled','b')
plot3(pyslicelst,pxlst,yTcheb,'r')
plot3(pyslicelst,pxlst,myfun(pxlst),'Black')
ylim([-1 1])
xlim([-1 1])

%{
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
 %}

else
%subplot(2,frm,k)
hold on 
scatter3(yslicelst,xlst,flst,'filled','b')
plot3(pyslicelst,pxlst,pylst,'r') 
plot3(pyslicelst,pxlst,myfun(pxlst),'Black')
plot3(pxlst,pyslicelst,pylst,'g') 
plot3(pxlst,pyslicelst,myfun(pxlst),'Black')
ylim([-1 1])
xlim([-1 1])

%{
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
%}
end

end
[X,Y] = meshgrid(a:2/24:b,a:2/24:b);
Z =my2dfun(X,Y);
surf(X,Y,Z);