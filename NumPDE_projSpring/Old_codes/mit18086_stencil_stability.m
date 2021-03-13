function mit18086_stencil_stability(x,nderiv,r)
%MIT18086_STENCIL_STABILITY(X,NDERIV,R)
%    Computes the stencil weights which approximate the
%    NDERIV-th derivative for a set of points given by X.
%    Also plots the von Neumann stability function of an
%    explicit time step method (with Courant number r)
%    solving the initial value problem u_t = u_nx, where
%    NDERIV space derivatives are taken.

% 02/2007 by Benjamin Seibold
% Feel free to modify for teaching and learning.

if nargin<3, r = 1; end
if nargin<2, nderiv = 3; end
if nargin<1, x = 0:3; end
V = rot90(vander(x));                      % Vandermonde matrix
b = zeros(size(x))'; b(nderiv+1) = factorial(nderiv);
s = V\b;                                      % stencil weights
disp('point positions:')
disp(sprintf('%7.2f',x))
disp(sprintf('stencil weights for %d. derivative:',nderiv))
disp(sprintf('%7.2f',s))
t = linspace(0,2*pi,200)';
g = 1+r*exp(i*t*x)*s;               % von Neumann growth factor
clf
patch(cos(t),sin(t),.9*[1 1 1])
hold on
plot(real(g),imag(g),'r-','linewidth',2)
hold off
axis equal
title('growth factor')