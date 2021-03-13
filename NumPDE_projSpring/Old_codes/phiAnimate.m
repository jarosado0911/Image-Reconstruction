function phiAnimate()
%Define domain
domain = [0 1 0 1];
nX = 2000;
nY = nX;

hx = (domain(2)-domain(1))/(nX-1);
hy = (domain(4)-domain(3))/(nY-1);

x = [0:hx:1];
y = [0:hy:1];

[X,Y]=meshgrid(x,y);
phi0 = zeros(nX,nY);


alpha=flip([0.1:0.1:1]);

figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.75]);
for k=1:length(alpha)
    snakeX=[];
    snakeY=[];
for i=1:nX
    for j=1:nY
        if (norm([X(i,j) Y(i,j)]-[0.5 0.5],2)<=0.5*alpha(k) && norm([X(i,j) Y(i,j)]-[0.5 0.5],2)>0.499*alpha(k))
            phi0(i,j)=1;
            snakeX=[snakeX X(i,j)];
            snakeY=[snakeY Y(i,j)];
        end
    end
end

scatter(snakeX,snakeY,'.b')
title(sprintf('\alpha = %f',alpha(k)));
xlabel('x')
ylabel('y')
xlim([-0.25 1.25]);
ylim([-0.25 1.25]);
legend('\phi_0(x,y)')
drawnow
end
end

