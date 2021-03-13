function [ret1,ret2] = pdcub(n,xyrange,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdcub".
%
% PadL = pdcub(n,xyrange)
% [cubature,PadL] = pdcub(n,xyrange,funct,opt1,opt2,...)
%
% This function computes the Padua points defined in the rectangle 
% [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)] and the 
% corresponding cubature weights as a matrix of abscissas (first
% column), ordinates (second column) and weights (third column).
% Optionally, it computes the cubature of the function funct (optional 
% input argument).
%-------------------------------------------------------------------------------
% INPUT.
%
% n        : interpolation degree
% xyrange  : a vector [a,b,c,d] defining the rectangle [a,b] x [c,d]
% funct    : function to be interpolated in the form 
%            funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%            optional arguments for funct 
%
% OUTPUT.
%
% PadL     : a matrix with the abscissas of Padua points in the 
%            first column, the ordinates in the second and the 
%            cubature weights in the third
% cubature : cubature of the integrand funct
%-------------------------------------------------------------------------------
% EXAMPLES.
%
% 1)
% Compute the Padua points and the cubature weights of degree 50
%
% xyrange = [0,1,0,1];
% PadL = pdcub(50,xyrange);
%
% 2)
% Compute the cubature of the Franke's function defined in funct.m
%
% xyrange = [0,1,0,1];
% [cubature,PadL] = pdcub(50,xyrange,@testfunct);
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% 1. pdpts
% 2. pdwtsMM
%-------------------------------------------------------------------------------
% Authors:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@math.unipd.it>   
%          Alvise Sommariva  <alvise@math.unipd.it>
%          Marco Vianello    <marcov@math.unipd.it>   
%
% Date: November 2009.
%-------------------------------------------------------------------------------

% Compute the cubature weights
PadL(:,3) = pdwtsMM(n,xyrange);
% Compute the Padua points in the rectangle defined by xyrange
PadL(:,1:2) = pdpts(n,xyrange);
if (nargin < 2)
  error('Too few input arguments')
elseif (nargin >= 3)
  funct = varargin{1};
  ret1 = PadL(:,3)'*feval(funct,PadL(:,1),PadL(:,2),varargin{2:end});
  ret2 = PadL;
else
% nargin == 2
  ret1 = PadL;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdcub
%
% at the Octave prompt
%
%!test
%! disp('Degree 0 cubature')
%! xyrange = [0,1,0,1];
%! [cubature,PadL] = pdcub(0,xyrange,@testfunct,7);
%! expected = 1;
%! assert(cubature/expected,1,10*eps);
%!test
%! disp('Franke''s function cubature')
%! xyrange = [0,1,0,1];
%! [cubature,PadL] = pdcub(500,xyrange,@testfunct);
%! expected = 4.069695894915573e-01;
%! assert(cubature/expected,1,50*eps);  
