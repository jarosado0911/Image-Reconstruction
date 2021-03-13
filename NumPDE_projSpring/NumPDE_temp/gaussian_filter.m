function E = gaussian_filter(raster,X,Y,r)

% Parameters
w=abs(X(1,2)-X(1,1)); %spacing of grid in the convolution
r=(r/w)*w;
%R = Matrix of correlation coefficents 
%W = weight kernel
E = raster; % input raster used instead of random generated white noise error. 

% correlation coefficients matrix R
m=(2*((2*r)/w))+1; % dimension of matrix R

% Gaussian autocorrelation function
% p(h)=exp(-3*(h^2/r^2)); %p(h) is the correlation coefficient of the DEM
% error at lag h and range r.
R = euklid_dist(m,m,w,w); % get euklid distance weight matrix & apply cellspacing w
R = exp( -3*( (R.^2) ./ (r^2) ) ); % apply corr coeff  as in paper

C = real(sqrtm(R.^2)); % this gives complex numbers, but close to 0, so ignore

s=sum(C(:).^2);
s = 1/sqrt(s);

W=s.*C;

E=conv2(E,W,'same');

function W = euklid_dist(msx,msy,wx,wy,n)
    msx=floor(msx/2);
    msy=floor(msy/2);
    [X,Y] = meshgrid(-msx:msx,-msy:msy);
    X=X*wx;
    Y=Y*wy;
    W = sqrt(X.^2+Y.^2);
%--------------------------------------------------------------------------
