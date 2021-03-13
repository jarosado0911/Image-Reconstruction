%-------------------------------------------------------------------------------
% This demo compares the relative errors (w.r.t. maximum deviation from average)
% of the interpolant of "f" in Padua Points computed respectively by FFT and 
% MM method.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdcfsFFT
% 2. pdcfsMM
% 3. pdpts
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

m=110;

%-------------------------- SETTINGS END HERE ----------------------------------

test_points=pdpts(m,xyrange);
f_test_points=feval(@testfunct,test_points(:,1),test_points(:,2));
aver_fvalue=mean(f_test_points);
max_dev_aver_fvalue=norm(f_test_points-aver_fvalue,inf);
    
for index=1:length(n_vett)
    
    n=n_vett(index);
    
    fprintf('\n \n \t [DEGREE]: %3.0f',n);
    
    % FFT METHOD.
    
    Pad = pdpts(n,xyrange);
    C0f_FFT = pdcfsFFT(Pad,@testfunct,xyrange);
    LnfX_FFT = pdval(C0f_FFT,xyrange,test_points);
  
    fft_err_inf=norm(f_test_points-LnfX_FFT,inf);
    fft_err_2rel(index)=fft_err_inf/max_dev_aver_fvalue;
    fprintf('\n \t [FFT] [REL.INTP.ERR.INF.]: %2.5e',fft_err_2rel(index));
    
    % MM METHOD.
    
    [X1,Y1,X2,Y2] = pdpts(n,xyrange);
    C0f_MM = pdcfsMM(X1,Y1,X2,Y2,@testfunct,xyrange);
    LnfX_MM = pdval(C0f_MM,xyrange,test_points);
    
    MM_err_inf=norm(f_test_points-LnfX_MM,inf);
    MM_err_2rel(index)=MM_err_inf/max_dev_aver_fvalue;
    
    fprintf('\n \t [MM ] [REL.INTP.ERR.INF.]: %2.5e',MM_err_2rel(index));
    
end

semilogy(n_vett,fft_err_2rel,'g-d',n_vett,MM_err_2rel,'b-x');
legend('FFT','MM');
xlabel('Degree');
ylabel('Interpolation relative error');

fprintf('\n');

