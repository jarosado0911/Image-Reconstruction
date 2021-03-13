clf
close all
clc
clear
my2dfun=@(x,y) (1+25.*x.^2.*y.^2).^(-1);
%my2dfun=@(x,y) sin(x).^2.*cos(y).^2;
a=-1;
b=1;

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

typePart=input('Type of points? [1] Equidistant or [2] Tchebyshev: ');
res=input('Resolution? ');
for j=1:res
clf
nxSamples=j+1;
nySamples=j+1;
hx=abs(b-a)/nxSamples;
hy=abs(b-a)/nySamples;

if typePart ~= 1
%Use Tchebyshev points
    yslice=(1/2)*(a+b)+(1/2)*(b-a)*cos( (2* (0:nySamples-1)+1) / (2*nySamples) * pi );
    xslice=(1/2)*(a+b)+(1/2)*(b-a)*cos( (2* (0:nxSamples-1)+1) / (2*nxSamples) * pi );
else
    yslice=linspace(a,b,nySamples);
    xslice=linspace(a,b,nxSamples);
end

%x_i sample points
xlst=xslice;

for k=1:length(yslice)
    %Fix y to get a slice
    myfun=@(x) my2dfun(x,yslice(k));
    
    %Get Sample f(x_i)=f_i
    flst=myfun(xlst);
    
    %Get newton coefficients
    clst=newtonCoeff(xlst,flst);
    
    %Get points for newton polynomial
    pxlst=linspace(a,b,1000);
    pylst=newtonP(pxlst,xlst,clst);
    
    yslicelst=yslice(k)*ones(1,length(xlst));
    pyslicelst=yslice(k)*ones(1,length(pxlst));
    
    hold on 
    %scatter3(xlst,yslicelst,flst,'filled','b')
    plot3(pxlst,pyslicelst,pylst,'r') 
    view(-30,50)
    %Plot of slice function
    %plot3(pxlst,pyslicelst,myfun(pxlst),'Black')
   
    
    ylim([a b])
    xlim([a b])
    zlim manual
    zlim([0 1])
end

%x_i sample points
ylst=yslice;

for k=1:length(xslice)
    %Fix y to get a slice
    myfun=@(y) my2dfun(xslice(k),y);
    
    %Get Sample f(x_i)=f_i
    flst=myfun(ylst);
    
    %Get newton coefficients
    clst=newtonCoeff(ylst,flst);
    
    %Get points for newton polynomial
    pylst=linspace(a,b,1000);
    pxlst=newtonP(pylst,ylst,clst);
    
    xslicelst=xslice(k)*ones(1,length(ylst));
    pxslicelst=xslice(k)*ones(1,length(pylst));
    
    %scatter3(xslicelst,ylst,flst,'filled','b')
    plot3(pxslicelst,pylst,pxlst,'r')
    view(-30,50)
    %Plot of slice function
    %plot3(pxslicelst,pylst,myfun(pylst),'Black')
    
    ylim([a b])
    xlim([a b])
    zlim manual
    zlim([0 1])
    
    
end
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
%Actual Plot
figure
[X,Y] = meshgrid(a:hx:b,a:hy:b);
Z =my2dfun(X,Y);
s=surf(X,Y,Z);
view(-30,50)