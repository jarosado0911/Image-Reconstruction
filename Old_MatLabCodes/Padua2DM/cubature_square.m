function [xv,yv,wv]=cubature_square(N,a,b,c,d,cubature_type)

%-------------------------------------------------------------------------------
% USAGE of "cubature_square".
%
% [x,y,w] = cubature_square(N,a,b,c,d,cubature_type)
%
% Compute the cubature nodes and weights of different tensor rules.
%
% INPUT.
% N             : number of nodes of 1D rule (and not the degree).
% a,b,c,d       : rectangle [a,b] x [c,d].
% cubature_type : 1 Fejer type 1 
%                 2 Fejer type 2 
%                 3 Clenshaw-Curtis 
%                 4 Gauss-Legendre 
%                 5 Gauss-Legendre-Lobatto.
%
% OUTPUT.
% x             : cubature nodes (abscissas)
% y             : cubature nodes (ordinates)
% w             : cubature weights
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Authors:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@math.unipd.it>   
%          Alvise Sommariva  <alvise@math.unipd.it>
%          Marco Vianello    <marcov@math.unipd.it>   
%
% Date: November 2009.
%-------------------------------------------------------------------------------

[nodes,weights]=quadrature_rules(N,cubature_type);

xnodes=(b+a)/2+(b-a)*nodes/2;
ynodes=(c+d)/2+(d-c)*nodes/2;

xweights=0.5*(b-a)*weights;
yweights=0.5*(d-c)*weights;

[x_mesh,y_mesh]=meshgrid(xnodes,xnodes);
[wx_mesh,wy_mesh]=meshgrid(xweights,yweights);
w_mesh=wx_mesh.*wy_mesh;

xv=x_mesh(:);
yv=y_mesh(:);
wv=w_mesh(:);


%-------------------------------------------------------------------------------
% 1. "quadrature_rules"
%-------------------------------------------------------------------------------
function [x,w]=quadrature_rules(n,quadrature_rule)

% NODES AND WEIGHTS OF SOME COMMON RULES, IN THE INTERVAL [-1,1].
% SEE WALDVOGEL PAPER.

switch quadrature_rule
    
case 1 % FEJER 1.
    [x,w]=fejer1(n); 
    
case 2 % FEJER 2.
    [x,w]=fejer2(n); 
    
case 3 % CLENSHAW-CURTIS.
    [x,w]=clenshaw_curtis(n);
    
case 4 % GAUSS-LEGENDRE.
    ab=r_jacobi(n,0,0);
    xw=gauss(n,ab);
    x=xw(:,1);
    w=xw(:,2);
    
case 5 % GAUSS-LEGENDRE-LOBATTO.
    ab=r_jacobi(n+2,0,0);
    xw=lobatto(n-2,ab,-1,1);
    x=xw(:,1);
    w=xw(:,2);    
    
end









%-------------------------------------------------------------------------------
% 2. fejer1
%-------------------------------------------------------------------------------

function [x,w]=fejer1(n)

%-------------------------------------------------------------------------------
% Waldvogel says that the nodes are 
%
%     x(k)=cos(theta(k)), k=1/2,\ldots,(n-1)/2
%
% where
%
%     theta(k)=k*pi/n,    k=1/2,\ldots,(n-1)/2
% 
% omitting the nodes x(0) and x(n).
%-------------------------------------------------------------------------------

% COMPUTING n-1 NODES.
k=(1:n)'-1/2;
x=cos(k*pi/n);


% Fejer1 nodes: k=1/2,3/2,...,n-1/2; vector of weights: wf1
N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';

v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
v1=v0(1:end-1)+conj(v0(end:-1:2)); 
w=ifft(v1);



%-------------------------------------------------------------------------------
% 3. fejer2
%-------------------------------------------------------------------------------

function [x,w]=fejer2(n)

%-------------------------------------------------------------------------------
% Waldvogel says that the nodes are 
%
%     x(k)=cos(theta(k)), k=0,\ldots,n
%
% where
%
%     theta(k)=k*pi/n,    k=0,\ldots,n
%-------------------------------------------------------------------------------
% IMPORTANT.
% The rule has n-1 points (index k=0 and k=n
% are dropped). 
% Consequently if one wants "n" points, 
% one has to set "n=n+1".
%-------------------------------------------------------------------------------

n=n+1;

% Computing weights.
k=(1:n-1)'; 
x=cos(k*pi/n);

% Computing nodes.
N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
v2=-v0(1:end-1)-v0(end:-1:2); 
w=ifft(v2);
w=w(2:end,1);



%-------------------------------------------------------------------------------
% 4. Clenshaw-Curtis
%-------------------------------------------------------------------------------

function [x,w]=clenshaw_curtis(n)

% Von Winkel code.

N1=n;
N=N1-1; bma=2;
c=zeros(N1,2);
c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
w=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
x=0.5*(N*bma*f(1:N1,2));



%-------------------------------------------------------------------------------
% 5. Gauss
%-------------------------------------------------------------------------------
%
% GAUSS Gauss quadrature rule.
%
%    Given a weight function w encoded by the nx2 array ab of the 
%    first n recurrence coefficients for the associated orthogonal
%    polynomials, the first column of ab containing the n alpha-
%    coefficients and the second column the n beta-coefficients, 
%    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
%    the n-point Gauss quadrature rule for the weight function w.
%    The nodes, in increasing order, are stored in the first 
%    column, the n corresponding weights in the second column, of
%    the nx2 array xw.
%
%    Supplied by Dirk Laurie; edited by Walter
%    Gautschi.
%
%-------------------------------------------------------------------------------
function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];



%-------------------------------------------------------------------------------
% 6. r_jacobi
%-------------------------------------------------------------------------------
%
% R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%
%    ab=R_JACOBI(n,a,b) generates the first n recurrence 
%    coefficients for monic Jacobi polynomials with parameters 
%    a and b. These are orthogonal on [-1,1] relative to the
%    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
%    are stored in the first column, the n beta-coefficients in
%    the second column, of the nx2 array ab. The call ab=
%    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
%    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
%
%    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%    Gautschi, 4-4-2002.
%
%-------------------------------------------------------------------------------

function ab=r_jacobi(N,a,b)
if nargin<2, a=0; end;  if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1, ab=[nu mu]; return, end 
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];



%-------------------------------------------------------------------------------
% 7. lobatto
%-------------------------------------------------------------------------------
%
% LOBATTO Gauss-Lobatto quadrature rule.
%
%    Given a weight function w encoded by the (n+2)x2 array ab of
%    the first n+2 recurrence coefficients for the associated
%    orthogonal polynomials, the first column of ab containing the
%    n+2 alpha-coefficients and the second column the n+2 beta-
%    coefficients, the call xw=LOBATTO(n,ab,endl,endr) generates 
%    the nodes and weights xw of the (n+2)-point Gauss-Lobatto 
%    quadrature rule for the weight function w having two
%    prescribed nodes endl, endr (typically the left and right end
%    points of the support interval of w, or points to the left
%    resp. to the right therof). The n+2 nodes, in increasing 
%    order, are stored in the first column, the n+2 corresponding
%    weights in the second column, of the (n+2)x2 array xw.
%        
%    For Jacobi weight functions, see also LOBATTO_JACOBI.
%    Supplied by Dirk Laurie; edited by Walter Gautschi.
%
%-------------------------------------------------------------------------------
function xw=lobatto(N,ab,endl,endr)
N0=size(ab,1); if N0<N+2, error('input array ab too short'), end
ab0=ab;
p0l=0; p0r=0; p1l=1; p1r=1;
for n=1:N+1
  pm1l=p0l; p0l=p1l; pm1r=p0r; p0r=p1r;
  p1l=(endl-ab0(n,1))*p0l-ab0(n,2)*pm1l;
  p1r=(endr-ab0(n,1))*p0r-ab0(n,2)*pm1r;
end
det=p1l*p0r-p1r*p0l;
ab0(N+2,1)=(endl*p1l*p0r-endr*p1r*p0l)/det;
ab0(N+2,2)=(endr-endl)*p1l*p1r/det;
xw=gauss(N+2,ab0);

