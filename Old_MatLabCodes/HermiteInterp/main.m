clf
clear
format long
myfun=@(x) (1+25*x.^2).^(-1);
mydfun=@(x) -(1+25*x.^2).^(-2).*(50.*x);
myddfun=@(x) 2*(1+25*x.^2).^(-3).*(50.*x).*(50.*x)-50*(1+25*x.^2).^(-2);
a=-1;
b=1;
nsamples=[4,6,9,11,14,16,19,21];
N=5;

xlst=linspace(a,b,N);
ylst=myfun(xlst);
ydlst=mydfun(xlst);
yddlst=myddfun(xlst);

%Number Of points
nPts=length(xlst);
%Order of highest derivative
mOrd=2;

%Degree of polynomial
polDeg=nPts*(mOrd+1)-1;

%Make a matrix of zeros
Amat=zeros(polDeg+1);

for n=1:mOrd+1
for k=1:nPts
    for j=n:polDeg+1
        Amat(k+nPts*(n-1),j)=xlst(k)^(j-n)*factorial(j-1)/factorial(j-1-n+1);
    end
end
end
Y=[ylst ydlst yddlst]';

VanMat=[Amat Y];
rowR=rref(VanMat);
coeffs=rowR(:,polDeg+2)';

xvals=linspace(a,b,1000);
yvals=polyNum(xvals,coeffs);

hold on
fplot(myfun,[a,b],'Black')
scatter(xlst,ylst,'r')
plot(xvals,yvals,'b')

%{
Bmat=zeros(polDeg+1);
for n=1:mOrd+1
for k=1:nPts
    for j=n:polDeg+1
        Bmat(k+nPts*(n-1),j)=factorial(j-1)/factorial(j-1-n+1);
    end
end
end
Bmat 
%}