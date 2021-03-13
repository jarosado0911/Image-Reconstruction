%Adaptive Integration Using Simpson
format long
clc
clear
clf
syms x
fun=@(x)(1+25.*x.^2).^(-1);
df1 = matlabFunction( diff(fun(x)) );
df2 = matlabFunction( diff(diff(fun(x))) );
df3 = matlabFunction( diff(diff(diff(fun(x)))) );
df4 = matlabFunction( diff(diff(diff(diff(fun(x))))));
a=-1;
b=1;
%tol=input('Enter a tolerance: ');
%tol=input('Enter a tolerance: ');

tol=[0.1,0.001,0.0001,0.00001,0.000001,0.000001];
frm=ceil(length(tol)/2);
for j=1:length(tol)
[myarea,plst,theError,theExactError,fevals]=adaptsimp(fun,a,b,tol(j),df4);
Iexact=integral(@(x)abs(fun(x)),a,b);
flst=fun(plst);
diff=abs(Iexact-myarea);
vec=zeros(length(plst),1);

subplot(2,frm,j)
hold on
fplot(fun,[a,b])
stem(plst,flst,'filled','r')
scatter(plst,vec,'filled','b')
ylim([0 1.25])
fprintf('For a tolerance of %f the approx area = %f.\n The exact area = %f.\n The approximate error = %d.\n The exact error(comp) = %d.\n The difference in integrals is %d.\n The number of function evaluations is %u.\n The number of points is %u.\n\n',tol(j), myarea,Iexact,theError,theExactError,diff,fevals,length(plst))
end
