function [trapArea,OutTrueError,ErrorEst] = TrapezoidArea(x1,x2,y1,y2,fun,df2)
trapArea=(1/2)*abs(x2-x1)*(abs(y1)+abs(y2));
OutTrueError=abs(integral(@(x)fun(x),x1,x2)-trapArea);
ErrorEst=abs((x2-x1)^3/12*fminbnd(@(x) -1*df2(x),x1,x2));
end

