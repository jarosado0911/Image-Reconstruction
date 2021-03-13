function E = gaussian_filter(raster,X,Y,r)
% -------------------------------------------------------------------
% USAGE: E = gaussian_filter(raster,X,Y,r)
% where: [raster] is the raster to apply filter to
%        [X] and [Y] are the coordinate vectors of [raster]
%        [r] is the (correlation) range of the error in units defined by 
%            the cellsize in X/Y
%
% returns matrix [E] of size [raster] = X * Y
% -------------------------------------------------------------------
% EXAMPLE: E = gaussian_filter(randn(100,100),[1:100],[1:100],5);

if nargin < 4
    error('This function needs 4 input parameters (inputraster,X,Y,range)');
end
% Parameters
w=abs(X(1,2)-X(1,1)); %spacing of grid in the convolution (resolution of DEM, as in the paper)
if w~=abs(Y(1,2)-Y(1,1))
    error('This approach only works for equidistant matrices (same cellsize in X and Y direction)!')
end
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
