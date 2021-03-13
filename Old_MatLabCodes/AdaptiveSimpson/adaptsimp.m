function [area,ptSet,errorTotal,error,fevals] = adaptsimp(fun,xi,xf,tol,df4)
[area,error,errorTotal]=simpson(xi,xf,fun,df4);
ptSet=[xi xf];
fevals=2;
if error>tol
    m=(xi+xf)/2;
    [a1, S1,E1,Ex1,feval1]=adaptsimp(fun,xi, m, tol,df4);
    [a2, S2,E2,Ex2,feval2]=adaptsimp(fun,m,xf, tol,df4);
    area=a1+a2;
    errorTotal=E1+E2;
    error=Ex1+Ex2;
    fevals=feval1+feval2;
    ptSet=[ptSet S1 S2];
    ptSet=sort(ptSet);
    ptSet=unique(ptSet','rows').';
end
end

