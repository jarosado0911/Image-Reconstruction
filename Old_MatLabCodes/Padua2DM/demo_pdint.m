%-------------------------------------------------------------------------------
% A demo that computes the values of the interpolant of a function "f" 
% in a prescribed set of points, knowing its samples at Padua Points.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdint
% 2. pdpts
% 3. funct
%-------------------------------------------------------------------------------
% TESTED ON:
%
% 1. GNU Octave 3.1.50.
% 2. Matlab 6.1.0.450. Release 12.1. 
%-------------------------------------------------------------------------------
% Authors: 
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@math.unipd.it>   
%          Alvise Sommariva  <alvise@math.unipd.it>
%          Marco Vianello    <marcov@math.unipd.it>   
%
% Date: November 2009.
%-------------------------------------------------------------------------------

clear all; close all; more off;

n=100;                      % PADUA POINTS DEGREE.
xyrange=[0,pi,-pi,pi];      % DEFINITION OF THE RECTANGLE.

M=30;                       % EVALUATION POINTS IN THE SQUARE: DEGREE OF PADUA
                            % POINTS USED FOR TESTS. 

%-------------------------- SETTINGS END HERE ----------------------------------

X = pdpts(M,xyrange);       % [xyrange(1),xyrange(2)]x[xyrange(3),xyrange(4)]
                            % STORED IN A MATRIX HAVING SIZE  "M x 2".
                            % THE FIRST COLUMN CONTAINS THE ABSCISSAS,
                            % THE SECOND COLUMN CONTAINS THE ORDINATES.
                                        
LnfX = pdint(n,xyrange,@testfunct,X); 
                            % "LnfX" IS THE INTERPOLANT AT PADUA POINTS OF 
                            % THE FUNCTION, STORED IN A MATRIX HAVING SIZE 
                            % "M x 1", "M" BEING "(n+1)*(n+2)/2".
                                        
fX=feval(@testfunct,X(:,1),X(:,2)); 
                            % EVALUATION OF THE FUNCTION AT PADUA POINTS.

absmaxerr=norm(fX-LnfX);    % COMPUTING MAXIMUM ABSOLUTE ERROR.

aver_fvalue=mean(fX);

index=find(abs(fX)>0);      % COMPUTING MAXIMUM RELATIVE ERROR.
relmaxerr=norm(fX(index)-LnfX(index),inf)/aver_fvalue;        

fprintf('\n\t [n]: %5.0f [max.abs.err.]: %2.2e [rel.abs.err.]: %2.2e \n \n',...
n,absmaxerr,relmaxerr);   

