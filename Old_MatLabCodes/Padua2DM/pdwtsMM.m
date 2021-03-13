function [ret1,L2,L] = pdwtsMM(n,xyrange)

%-------------------------------------------------------------------------------
% USAGE of "pdwtsMM".
%
% L = pdwtsMM(n)
% L = pdwtsMM(n,xyrange)
% [L1,L2] = pdwtsMM(n)
% [L1,L2] = pdwtsMM(n,xyrange)
% [L1,L2,L] = pdwtsMM(n)
% [L1,L2,L] = pdwtsMM(n,xyrange)
%
% Compute the cubature weights L by MM so that, if Pad is the matrix of Padua
% points computed through the call Pad = pdpts(n) or Pad = pdpts(n,xyrange), 
% then the cubature of the function funct is given by 
%
% L'*funct(Pad(:,1),Pad(:,2))
%
% An alternative function for the computation of the cubature weights
% by FFT is pdwtsFFT.
%
% Otherwise, one can compute the cubature weights L1 and L2, so that, 
% if X1,Y1 and X2,Y2 are the subgrids of Padua points computed through 
% the call [X1,Y1,X2,Y2] = pdpts(n) or [X1,Y1,X2,Y2] = pdpts(n,xyrange),
% then the cubature of the function funct is given by 
%
% sum(sum(L1.*funct(X1,Y1)))+sum(sum(L2.*funct(X2,Y2)))
%-------------------------------------------------------------------------------
% INPUT.    
%
% n       : interpolation degree
% xyrange : an optional vector [a,b,c,d] defining the rectangle 
%           [a,b] x [c,d]. By default, xyrange = [-1,1,-1,1]
%
% OUTPUT.   
%
% L       : cubature weights associated to the matrix Pad of Padua points
% L1,L2   : cubature weights associated to the subgrids 
%           X1,Y1 and X2,Y2 of Padua points
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% no external function is used.
%-------------------------------------------------------------------------------
% Authors:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@math.unipd.it>   
%          Alvise Sommariva  <alvise@math.unipd.it>
%          Marco Vianello    <marcov@math.unipd.it>   
%
% Date: November 2009.
%-------------------------------------------------------------------------------

if (nargin == 1)
  xyrange = [-1,1,-1,1];
end
if (n == 0)
% degree 0
  L1 = (xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3));
  L2 = zeros(1,0);
  L = L1;
else
  argn1 = linspace(0,pi,n+1);
  argn2 = linspace(0,pi,n+2);
  k = [0:2:n]';
  l = (n-mod(n,2))/2+1;
% even-degree Chebyshev polynomials on the subgrids
  TE1 = cos(k*argn1(1:2:n+1));
  TE1(2:l,:) = TE1(2:l,:)*sqrt(2);
  TO1 = cos(k*argn1(2:2:n+1));
  TO1(2:l,:) = TO1(2:l,:)*sqrt(2);
  TE2 = cos(k*argn2(1:2:n+2));
  TE2(2:l,:) = TE2(2:l,:)*sqrt(2);
  TO2 = cos(k*argn2(2:2:n+2));
  TO2(2:l,:) = TO2(2:l,:)*sqrt(2);
% even,even moments matrix
  mom = 2*sqrt(2)./(1-k.^2);
  mom(1) = 2;
  [M1,M2] = meshgrid(mom);
  M = M1.*M2;
  M0 = fliplr(triu(fliplr(M)));
% interpolation weights matrices
  W1 = 2*ones(l)/(n*(n+1));
  W2 = 2*ones((n+mod(n,2))/2+1,(n+mod(n,2))/2)/(n*(n+1));
  W1(:,1) = W1(:,1)/2;
  W2(1,:) = W2(1,:)/2;
  if (mod(n,2) == 0)
    M0(n/2+1,1) = M0(n/2+1,1)/2;
    W1(:,n/2+1) = W1(:,n/2+1)/2;
    W1(n/2+1,:) = W1(n/2+1,:)/2;
  else
    W2((n+1)/2+1,:) = W2((n+1)/2+1,:)/2;
    W2(:,(n+1)/2) = W2(:,(n+1)/2)/2;
  end
% cubature weights as matrices on the subgrids.
  L1 = W1.*(TE1'*M0*TO2)';
  L2 = W2.*(TO1'*M0*TE2)';
  L1 = L1*(xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))/4;
  L2 = L2*(xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))/4;
  if (nargout ~= 2)
    if (mod(n,2) == 0)
      L = zeros(n/2+1,n+1);
      L(:,1:2:n+1) = L1;
      L(:,2:2:n+1) = L2;
      L = L(:);
    else
      L = zeros((n+1)/2,(n+2));
      L = [L1',L2']';
      L = L(:);
    end
  end
end
if (nargout == 0 || nargout == 1)
  ret1 = L;
else
  ret1 = L1;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdwtsMM
%
% at the Octave prompt
%
%!test
%! disp('Degree 0 weight')
%! xyrange = [0,1,0,1];
%! L = pdwtsMM(0,xyrange);
%! expected = 1;
%! assert(L/expected,1,10*eps);
%!test
%! disp('Degree 0 weight, subgrids')
%! xyrange = [-1,1,-1,1];
%! [L{1},L{2}] = pdwtsMM(0,xyrange);
%! expected{1} = 4; 
%! expected{2} = zeros(1,0);
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 1 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsMM(1,xyrange);
%! expected = [0.5;0.25;0.25];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 1 weights, subgrids')
%! xyrange = [-1,1,-1,1];
%! [L{1},L{2}] = pdwtsMM(1,xyrange);
%! expected{1} = 2;
%! expected{2} = [1;1];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 2 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsMM(2,xyrange);
%! expected = [1/6;0;1/9;5/9;1/6;0];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 2 weights, subgrids')
%! xyrange = [-1,1,-1,1];
%! [L{1},L{2}] = pdwtsMM(2,xyrange);
%! expected{1} = [2/3,2/3;0,0];
%! expected{2} = [4/9;20/9];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 1 weights, subgrids, xyrange')
%! xyrange = [0,1,0,1];
%! [L{1},L{2}] = pdwtsMM(1,xyrange);
%! expected{1} = 0.5;
%! expected{2} = [0.25;0.25];
%! assert(L,expected,10*eps);
