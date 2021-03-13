function [phi0,snakeX,snakeY,x,y]=initial_sqr(npts,r,grph)

%Define domain
domain = [0 1 0 1];
nX = npts;
nY = nX;
eps=0.96;
r1=r;
r2 = eps*r;

hx = (domain(2)-domain(1))/(nX-1);
hy = (domain(4)-domain(3))/(nY-1);

x = (0:hx:1);
y = (0:hy:1);

[X,Y]=meshgrid(x,y);
phi0 = zeros(nX,nY);

snakeX=[];
snakeY=[];

r = 40;

for i=1:nX
    for j=1:nY
        if i==40 || i==length(x)-40
            phi0(i,j)=1;
            snakeX=[snakeX X(i,j)];
            snakeY=[snakeY Y(i,j)];
        elseif j==40 || j==length(x)-40
            phi0(i,j)=1;
            snakeX=[snakeX X(i,j)];
            snakeY=[snakeY Y(i,j)];
        end
    end
end

if grph ==1
    figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.75]);
scatter(snakeX,snakeY,'.b')
xlim([-0.25 1.25]);
ylim([-0.25 1.25]);
title('Initial Circle Snake');
xlabel('x')
ylabel('y')

legend('\phi_0(x,y)')
end
end

