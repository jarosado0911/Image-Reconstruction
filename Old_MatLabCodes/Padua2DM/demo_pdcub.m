%-------------------------------------------------------------------------------
% A demo that shows how to compute the cubature of a function "f" on a rectangle
% defined by "xyrange".
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdcub
% 2. funct
% 3. pdwtsFFT
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

clear all; close all; 
%more off;

n=100;                      % PADUA POINTS DEGREE.
xyrange=[0,1,0,1];          % DEFINITION OF THE RECTANGLE.

f=inline(['.75*exp(-((9*x-2).^2 + (9*y-2).^2)/4)+',...
         '.75*exp(-((9*x+1).^2)/49 - (9*y+1)/10)+',...
         '.5*exp(-((9*x-7).^2 + (9*y-3).^2)/4)-',...
         '.2*exp(-(9*x-4).^2 - (9*y-7).^2)']);   

n_exact=200;                % DEGREE OF PADUA POINTS USED FOR COMPUTING AN
                            % APPROX. OF THE EXACT INTEGRAL.

%-------------------------- SETTINGS END HERE ----------------------------------

% COMPUTING APPROXIMATION OF EXACT INTEGRAL.
PadL_exact = pdpts(n_exact,xyrange);
fPadL_exact=feval(f,PadL_exact(:,1),PadL_exact(:,2));
w_exact = pdwtsFFT(n_exact,xyrange);
I_exact = w_exact'*fPadL_exact;   

% FIRST METHOD.
t(1)=cputime;
PadL = pdcub(n,xyrange);
X1=PadL(:,1);
X2=PadL(:,2);
W =PadL(:,3);

f_PadL=feval(f,X1,X2);
I=W'*f_PadL;
t(2)=cputime;   

fprintf('\n \t >> FIRST METHOD.'); 
fprintf('\n \t [DEGREE PADUA POINTS]: %5.0f',n); 
fprintf('\n \t [NUMBER PADUA POINTS]: %5.0f',length(W)); 
fprintf('\n \t [CUB.RES.PD.  ]: %20.20e',I); 
fprintf('\n \t [CUB.RES.     ]: %20.20e',I_exact); 
fprintf('\n \t [CPUTIME 1-ST METHOD]: %2.2e',t(2)-t(1)); 


% SECOND METHOD.

t(1)=cputime;
[cubature,PadL] = pdcub(n,xyrange,f);
t(2)=cputime;

fprintf('\n \n \t >> SECOND METHOD.');
fprintf('\n \t [DEGREE PADUA POINTS]: %5.0f',n); 
fprintf('\n \t [NUMBER PADUA POINTS]: %5.0f',size(PadL,1)); 
fprintf('\n \t [CUB.RES.PD.  ]: %20.20e',cubature); 
fprintf('\n \t [CUB.RES.     ]: %20.20e',I_exact); 
fprintf('\n \t [CPUTIME 2-ND METHOD]: %2.2e',t(2)-t(1)); 

% FUNCTION SAVED ON A FILE "funct.m".

t(1)=cputime;
[cubature,PadL] = pdcub(n,xyrange,@testfunct);
t(2)=cputime;

fprintf('\n \n \t >> THIRD METHOD.');
fprintf('\n \t [DEGREE PADUA POINTS]: %5.0f',n); 
fprintf('\n \t [NUMBER PADUA POINTS]: %5.0f',size(PadL,1)); 
fprintf('\n \t [CUB.RES.PD.  ]: %20.20e',cubature); 
fprintf('\n \t [CUB.RES.     ]: %20.20e',I_exact); 
fprintf('\n \t [CPUTIME 3-RD METHOD]: %2.2e',t(2)-t(1)); 

% FUNCTION SAVED ON A FILE "testfunct.m" WITH ADDITIONAL PARAMETERS 
% (e.g. TO CHOOSE THE FUNCTION IN A SET OF FUNCTIONS, SEE "testfunct.m").

t(1)=cputime;
[cubature,PadL] = pdcub(n,xyrange,@testfunct,1);
t(2)=cputime;

fprintf('\n \n \t >> FOURTH METHOD.');
fprintf('\n \t [DEGREE PADUA POINTS]: %5.0f',n); 
fprintf('\n \t [NUMBER PADUA POINTS]: %5.0f',size(PadL,1)); 
fprintf('\n \t [CUB.RES.PD.  ]: %20.20e',cubature); 
fprintf('\n \t [CUB.RES.     ]: %20.20e',I_exact); 
fprintf('\n \t [CPUTIME 4-TH METHOD]: %2.2e \n \n',t(2)-t(1)); 
