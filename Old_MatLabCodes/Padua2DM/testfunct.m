function z = testfunct(x,y,iindex)

%-------------------------------------------------------------------------------
% USAGE of "testfunct".
%
% z = testfunct(x,y)
% z = testfunct(x,y,iindex)
%
% Default. Evaluate the Franke's function or the bivariate function defined 
% by iindex. If varargin is choosen, a different function can be used.
%-------------------------------------------------------------------------------
% INPUT.
% x,y    : evaluate the function in the points (x,y)
% iindex : index of the various functions. By default, iindex = 1.
%
% OUTPUT.
% z      : value of the function in the points (x,y)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% calls    : no external routine is used by this code
%-------------------------------------------------------------------------------
% Authors:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@math.unipd.it>   
%          Alvise Sommariva  <alvise@math.unipd.it>
%          Marco Vianello    <marcov@math.unipd.it>   
%
% Date: November 2009.
%-------------------------------------------------------------------------------

if (nargin == 2)
  iindex = 1;
end

switch iindex
case 1
% 1. Franke function.
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 2.1547794245591083e+000 with an estimated absolute error of
% 8.88e-016.
%? 3. The value of the definite integral on the square [0,1] x [0,1], 
%? obtained using a Padua Points cubature formula of degree 500, 
%? is 4.06969589491556e-01 with an estimated absolute error of
%? 8.88e-016.
%  Maple: 0.40696958949155611906
  z = 3/4*exp(-((9*x-2).^2+(9*y-2).^2)/4) +...
      3/4*exp(-(9*x+1).^2/49-(9*y+1)/10) + ...
      1/2*exp(-((9*x-7).^2+(9*y-3).^2)/4) -...
      1/5*exp(-(9*x-4).^2-(9*y-7).^2);
case 2
% 1. The value of the definite integral on the square [-1,1] x [-1,1],
% obtained using a Padua Points cubature formula of degree 2000, 
% is 3.9129044444568244e+000 with an estimated absolute error of 3.22e-010.
  z = ((x-0.5).^2+(y-0.5).^2).^(1/2);
case 3
% 1. Bivariate polynomial having moderate degree.
% 2. The value of the definite integral on the square [-1,1] x
% [-1,1], up to machine precision, is 18157.16017316017 (see ref. 6).
% 3. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 1.8157160173160162e+004.
% 4. 2D modification of an example by L.N.Trefethen (see ref. 7), f(x)=x^20.
  z = (x+y).^20;                     
case 4
% 1. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 2000, 
% is 2.1234596326670683e+001 with an estimated absolute error of 7.11e-015.
  loc_arg = (x-0.5).^2+(y-0.5).^2;
  z = exp(loc_arg);
case 5
% 1. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 2000, 
% is 3.1415926535849605e-002 with an estimated absolute error of 3.47e-017.
  loc_arg = (x-0.5).^2+(y-0.5).^2;
  z = exp(-100*loc_arg);
case 6
% 1. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 4.3386955120336568e-003 with an estimated absolute error of 2.95e-017.
  z = cos(30*(x+y));
case 7
% 1. Constant. To test interpolation and cubature at degree 0.
% 2. The value of the definite integral on the square [-1,1] x [-1,1]
% is 4.  
  z = ones(size(x));
case 8 
% 1. The value of the definite integral on the square [-1,1] x [-1,1] 
% is up to machine precision is 5.524391382167263 (see ref. 6).
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 5.5243913821672628e+000 with an estimated absolute error of 0.00e+000.
% 3. 2D modification of an example by L.N.Trefethen (see ref. 7), 
% f(x)=exp(x).
  z = exp(x+y);                      
case 9
% 1. Bivariate Runge function: as 1D complex function is analytic 
% in a neighborhood of [-1; 1] but not throughout the complex plane.
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% up to machine precision, is 0.597388947274307 (see ref. 6).
% 3. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 5.9738894727430725e-001 with an estimated absolute error of 0.00e+000.
% 4. 2D modification of an example by L.N.Trefethen (see ref. 7), 
% f(x)=1/(1+16*x^2).
  den = 1+16*(x.^2+y.^2); 
  z = 1./den;  
case 10
% 1. Low regular function.
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% up to machine precision, is 2.508723139534059 (see ref. 6).
% 3. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 2.5087231395340579e+000 with an estimated absolute error of 0.00e+000.
% 4. 2D modification of an example by L.N.Trefethen (see ref. 7), 
% f(x)=abs(x)^3.
  z = (x.^2+y.^2).^(3/2);            
case 11
% 1. Bivariate gaussian: smooth function.
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% up to machine precision, is 2.230985141404135 (see ref. 6).
% 3. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 500, 
% is 2.2309851414041333e+000 with an estimated absolute error of 2.66e-015.
% 4. 2D modification of an example by L.N.Trefethen (see ref. 7), 
% f(x)=exp(-x^2).
  z = exp(-x.^2-y.^2);               
case 12
% 1. Bivariate example stemming from a 1D C-infinity function.
% 2. The value of the definite integral on the square [-1,1] x [-1,1], 
% up to machine precision, is 0.853358758654305 (see ref. 6).
% 3. The value of the definite integral on the square [-1,1] x [-1,1], 
% obtained using a Padua Points cubature formula of degree 2000, 
% is 8.5335875865430544e-001 with an estimated absolute error of 3.11e-015.
% 4. 2D modification of an example by L.N.Trefethen (see ref. 7), 
% f(x)=exp(-1/x^2).
  arg_z = (x.^2+y.^2);
% Avoid cases in which "arg_z=0", setting only in those instances 
% "arg_z=eps".
  arg_z = arg_z+(1-abs(sign(arg_z)))*eps; 
  arg_z = 1./arg_z;
  z = exp(-arg_z);         
otherwise
  error('invalid function index')
end


% Octave testing: type
%
% test testfunct
%
% at the Octave prompt
%
%!test
%!  
%! disp('Franke''s function')
%! expected = 7.664205912849231e-01;
%! assert(testfunct(0,0)/expected,1,10*eps)
