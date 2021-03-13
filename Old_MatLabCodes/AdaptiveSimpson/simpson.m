function [area,error,ErrorEst] = simpson(x1,x2,fun,df4)
h=abs(x2-x1)/2;
xm=(x2+x1)/2;
area=(h/3)*(fun(x1)+4*fun(xm)+fun(x2));
error=abs(integral(@(x)fun(x),x1,x2)-abs(area));
ErrorEst=abs((x2-x1)^5/90*fminbnd(@(x) -1*df4(x),x1,x2));
end

