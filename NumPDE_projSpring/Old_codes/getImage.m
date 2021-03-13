function getImage()
RGB = imread('grumpycat.png');

u_initial = double(rgb2gray(RGB)./200);
mySize = min(size(u_initial));
u_initial=u_initial(1:mySize,1:mySize);
size(u_initial)
x=linspace(0,1,mySize);
imagesc(x,x,u_initial)
colorbar
caxis([0 1])
end