function u = getImage()
RGB = imread('dend7.png');

u = double(rgb2gray(RGB)./200);
mySize = floor(min(size(u)));
u=u(1:mySize,1:mySize);
u=abs(u-1);
end