function [Dx,Dy] = grad(M,dx,dy)

[r,c]=size(M);
Dx = zeros(r,c);
Dy = zeros(r,c);

for i=2:r-1
    for j=2:c-1
        Dx(i,j)=(M(i,j+1)-M(i,j-1))/(2*dx);
        Dy(i,j)=(M(i+1,j)-M(i-1,j))/(2*dy);
    end
end

for i=2:r-1
    Dx(i,1)=(M(i,2)-M(i,1))/dx;
    Dx(i,end)=(-M(i,r-1)+M(i,end))/dx;
    
    Dy(1,i)=(-M(2,i)+M(1,i))/dy;
    Dy(end,i)=(M(r-1,i)-M(end,i))/dy;
end

end

