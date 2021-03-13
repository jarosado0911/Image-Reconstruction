function s=mit18086_stencil_stability(x,nderiv)
if nargin<2, nderiv = 3; end
if nargin<1, x = 0:3; end
V = rot90(vander(x));                      % Vandermonde matrix
b = zeros(size(x))'; b(nderiv+1) = factorial(nderiv);
s = V\b;                                      % stencil weights
