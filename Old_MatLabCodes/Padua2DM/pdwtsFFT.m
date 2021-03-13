function L = pdwtsFFT(n,xyrange)

%-------------------------------------------------------------------------------
% USAGE of "pdwtsFFT".
%
% L = pdwtsFFT(n)
% L = pdwtsFFT(n,xyrange)
%
% Compute the cubature weights L by FFT so that, if Pad is the matrix of Padua 
% points computed through the call Pad = pdpts(n) or Pad = pdpts(n,xyrange), 
% then the cubature of the function fun is given by 
%
% L'*fun(Pad(:,1),Pad(:,2))
%
% The recommended function for the computation of the cubature weights
% by matrix multiplications is pdwtsMM.
%
%-------------------------------------------------------------------------------
% INPUT.    
%
% n       : interpolation degree
% xyrange : an optional vector [a,b,c,d] defining the rectangle 
%           [a,b] x [c,d]. By default, xyrange = [-1,1,-1,1]
%
% OUTPUT.   
%
% L       : cubature weights
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% no external function is used.
%-------------------------------------------------------------------------------
% Author:  
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
  L = (xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3));
else
% even,even moments matrix
  k = [2:2:n]';
  mom = zeros(n+1,1);
  mom(1) = 2*sqrt(2);
  mom(3:2:n+1) = 4./(1-k.^2);
  [M1,M2] = meshgrid(mom);
  M = M1.*M2;
  M(1,n+1) = M(1,n+1)/2;
  M(1,:) = M(1,:)/sqrt(2);
  M(:,1) = M(:,1)/sqrt(2);
  M0 = fliplr(triu(fliplr(M)));
  M0hat = real(fft(M0,2*(n+1)));
%  M0hat = M0hat(1:n+2,:);
  M0hathat = real(fft(M0hat(1:n+2,:),2*n,2));
  TMT = M0hathat(:,1:n+1);
% interpolation weights
  W = 2*ones(n+2,n+1)/(n*(n+1));
  W(1,:) = W(1,:)/2;
  W(:,1) = W(:,1)/2;
  W(n+2,:) = W(n+2,:)/2;
  W(:,n+1) = W(:,n+1)/2;
  L = W.*TMT;
  findM = [2:2:(n+1)*(n+2)]';
  if (mod(n,2) == 0)
    add = repmat([0,1],(n+2)/2,n/2);
    add = [add,zeros((n+2)/2,1)];
    findM = findM-add(:);
  end
%  [M1,M2] = meshgrid([0:n],[0:n+1]);
%  findM = find(mod(M1+M2,2));
% cubature weights
  L = L(findM);
  L = L*(xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))/4;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdwtsFFT
%
% at the Octave prompt
%
%!test
%! disp('Degree 0 weight')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(0,xyrange);
%! expected = 1;
%! assert(L/expected,1,10*eps);
%!test
%! disp('Degree 1 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(1,xyrange);
%! expected = [0.5;0.25;0.25];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 2 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(2,xyrange);
%! expected = [1/6;0;1/9;5/9;1/6;0];
%! assert(L,expected,10*eps);
