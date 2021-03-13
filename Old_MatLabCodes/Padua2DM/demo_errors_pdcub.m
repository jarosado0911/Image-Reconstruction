%-------------------------------------------------------------------------------
% This demo computes the cubature of a function "f" on a rectangle by an 
% algebraic formula of degree "n" on Padua points. 
% A comparison of relative errors between FFT method and MM method is made.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdpts
% 2. pdwtsFFT
% 3. pdwtsMM
% 4. funct
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
% more off;

n_vett=2:2:100;                            % PADUA POINTS DEGREES VECTOR.
xyrange=[0,1,0,1];                         % DEFINITION OF THE RECTANGLE.
n_exact=200;

%-------------------------- SETTINGS END HERE ----------------------------------

PadL_exact = pdpts(n_exact,xyrange);
fPadL_exact=feval(@testfunct,PadL_exact(:,1),PadL_exact(:,2));
w_exact = pdwtsFFT(n_exact,xyrange);
I_exact = w_exact'*fPadL_exact;                                           

for index=1:length(n_vett)
    
    n=n_vett(index);
    
    PadL = pdpts(n,xyrange);
    fPadL=feval(@testfunct,PadL(:,1),PadL(:,2));
    
    % FFT METHOD.
    wFFT = pdwtsFFT(n,xyrange);
    cubFFT(index) = wFFT'*fPadL;
    relerrFFT(index)=abs((cubFFT(index)-I_exact)/I_exact);
    
    % MM METHOD.
    wMM = pdwtsMM(n,xyrange);
    cubMM(index) = wMM'*fPadL;
    relerrMM(index)=abs((cubMM(index)-I_exact)/I_exact);
    
 fprintf('\n \n \t [DEGREE]: %4.0f [REL.ERR.FFT]: %2.5e [REL.ERR.MM]: %2.5e\n',n,relerrFFT(index),relerrMM(index));
    
end

semilogy(n_vett,relerrFFT,'g-d',n_vett,relerrMM,'b-x');
legend('FFT','MM');
xlabel('Degree')
ylabel('Relative cubature error')

fprintf('\n \n');