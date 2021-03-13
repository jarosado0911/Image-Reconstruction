function [area,ptSet,errorTotal,error,fevals] = adaptInt(fun,xi, xf, tol,df2)
y1=fun(xi);
y2=fun(xf);
[area,error,errorTotal]=TrapezoidArea(xi,xf,y1,y2,fun,df2);
fevals=2;
ptSet=[xi xf];
if error>=tol
    m=(xi+xf)/2;
    [a1, S1,E1,Ex1,feval1]=adaptInt(fun,xi, m, tol,df2);
    [a2, S2,E2,Ex2,feval2]=adaptInt(fun,m,xf, tol,df2);
    area=a1+a2;
    errorTotal=E1+E2;
    error=Ex1+Ex2;
    fevals=feval1+feval2;
    ptSet=[ptSet S1 S2];
    ptSet=sort(ptSet);
    ptSet=unique(ptSet','rows').';
end
end

