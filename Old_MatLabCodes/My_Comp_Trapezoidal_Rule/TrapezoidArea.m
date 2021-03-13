function [trapArea,OutError,OutTrueError] = TrapezoidArea(x1,x2,y1,y2,dfun,fun)
trapArea=(1/2)*abs(x2-x1)*(abs(y1)+abs(y2));
OutError=(abs(x1-x2)^2/12)*(dfun(x1)-dfun(x2));
OutTrueError=integral(fun,x1,x2)-trapArea;
end

