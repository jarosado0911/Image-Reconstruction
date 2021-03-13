function [u_init,x,y]=u0(npts,noise)

%Define domain
domain = [0 1 0 1];
nX = npts;
nY = nX;

hx = (domain(2)-domain(1))/(nX-1);
hy = (domain(4)-domain(3))/(nY-1);

x = (0:hx:1);
y = (0:hy:1);

[X,Y]=meshgrid(x,y);
u_init = zeros(nX,nY);
for i=1:nX
    for j=1:nX
        if (norm([X(i,j) Y(i,j)]-[0.42 0.42],inf)<=0.1 || norm([X(i,j) Y(i,j)]-[0.55 0.55],2)<=0.1)
            u_init(i,j)=1;
        end
    end
end
    u_init = u_init + rand(npts,npts)*noise;
end

