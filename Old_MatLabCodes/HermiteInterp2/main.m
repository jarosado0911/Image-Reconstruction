%MAIN HERMITE INTERPOLATION
clf
close all
clear
clc
format long

myfun=@(x) (1+25*x.^2).^(-1);
mydfun=@(x) -(1+25*x.^2).^(-2).*(50.*x);
myddfun=@(x) 2*(1+25*x.^2).^(-3).*(50.*x).*(50.*x)-50*(1+25*x.^2).^(-2);
a=-1;
b=1;
nsamples=[1:1:20];
choose=input('Equidistant [1] or Tchebyshev[2]? ');
frm=ceil(length(nsamples)/2);

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

for k=1:length(nsamples)
    clf
N=nsamples(k);
if choose == 1
    xlst=linspace(a,b,N);
    ylst=myfun(xlst);
    ydlst=mydfun(xlst);
    yddlst=mydfun(xlst);
elseif choose == 2
    xlst=(1/2)*(a+b)+(1/2)*(b-a)*cos( (2* (0:N-1)+1) / (2*N) * pi );
    ylst=myfun(xlst);
    ydlst=mydfun(xlst);
    yddlst=mydfun(xlst);
end
%Number Of points
nPts=length(xlst);
%Order of highest derivative
mOrd=2;

%Degree of polynomial
polDeg=nPts*(mOrd+1)-1;

Xlst=[];
Ylst=[];
YDlst=[];
YDDlst=[];
for n=1:nPts
    Xlst=[Xlst,repmat(xlst(n),mOrd+1,1)'];
    Ylst=[Ylst,repmat(ylst(n),mOrd+1,1)'];
    YDlst=[YDlst,repmat(ydlst(n),mOrd+1,1)'];
    YDDlst=[YDDlst,repmat(yddlst(n)./2,mOrd+1,1)'];
end

coeffs = hermiteInterp(Xlst,Ylst,YDlst,YDDlst,polDeg);

xvals=linspace(-1,1,1000);
yvals=newtonP(xvals,Xlst(1:polDeg+1),coeffs);

%subplot(2,frm,k)
hold on
plot(xvals,yvals,'b')
scatter(xlst,ylst,'r')
fplot(myfun,[a,b],'Black')
%pbaspect([2 2 2])
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

%PROBLEM #4
%For equidistant points, there is no convergence at the endpoints for
%Hermite interpolation.

%For Tchebyshev points there is convergence, but I noticed some disparity
%at x=-1 for more points.
