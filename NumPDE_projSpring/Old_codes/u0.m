function [u_init,x,y]=u0(npts)

%Define domain
domain = [0 1 0 1];
nX = npts;
nY = nX;

hx = (domain(2)-domain(1))/(nX-1);
hy = (domain(4)-domain(3))/(nY-1);

x = (0:hx:1);
y = (0:hy:1);

[X,Y]=meshgrid(x,y);
u_init = ones(nX,nY);

for i=1:nX
    for j=1:nX
        if (norm([X(i,j) Y(i,j)]-[0.45 0.40],1)<=0.15 ||...
                norm([X(i,j) Y(i,j)]-[0.65 0.45],5)<=0.1||...
                norm([X(i,j) Y(i,j)]-[0.5 0.625],2)<=0.1)
            u_init(i,j)=0;
        end
        
%         if (norm([X(i,j) Y(i,j)]-[0.2 0.7],2)<=0.09 &&...
%                 norm([X(i,j) Y(i,j)]-[0.2 0.7],2)>=0.025)
%             u_init(i,j)=0.25;
%         end
    end
end

end

