format long
clf
infun = input('Give function in x: ','s');
indfun= input('Give the derivative in x: ','s');
inddfun= input('Give the second derivative in x: ','s');

%Number of point N
N=input('How many points? ','s');
N=str2double(N);

%The function and its first two derivatives
fun = str2func(['@(x)',infun]);
dfun= str2func(['@(x)',indfun]);
ddfun =str2func(['@(x)',inddfun]);

%Domain
a=input('Start at a = ','s');
b=input('End at b = ','s');
a=str2double(a);
b=str2double(b);

%Uniform partition width
h=abs(b-a)/N;

%Make a list of x-samples
xlst=linspace(a,b,N+1);

%Get f list, here I am using x^2 function
flst=fun(xlst);

%Generate Trapezoid Weights
wlst=ones(1,N+1);
wlst(1)=1/2;
wlst(N+1)=1/2;

%Calculate Integral using trapezoidal rule
Th=h*dot(abs(flst),wlst);

%MatLab Calculation
posfun=@(x) abs(fun(x));
Iexact=integral(posfun,a,b);

%Get max error and compensate using error
xEmax=fminbnd(@(x)-ddfun(x), 0, 1);
ErrTmax=h^2*(b-a)/12*ddfun(xEmax);
comp1=Th-ErrTmax;

%Compensate Integral using correction
ErrComp=h^2/12*(dfun(a)-dfun(b));
comp2=Th+ErrComp;

fplot(fun,[a,b],'b')
hold on
for k=1:1:length(xlst)-1
    rx = [xlst(k) xlst(k) xlst(k+1) xlst(k+1)];
    ry = [flst(k) 0 flst(k+1) 0];
    k = convhull(rx, ry);
    fill (rx(k), ry(k), 'g','facealpha', 0.23); 
end
stem(xlst,flst,'r','filled')
plot(xlst,flst,'r')