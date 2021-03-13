%-------------------------------------------------------------------------------
% This demo computes the average CPU time needed to compute the cubature weights
% relative to an algebraic formula of degree "n" on Padua points. 
% A comparison between FFT method and MM method is made.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdwtsFFT
% 2. pdwtsMM
%-------------------------------------------------------------------------------
% TESTED ON:
%
% 1. GNU Octave 3.1.50.
% 2. Matlab 6.1.0.450. Release 12.1.
%-------------------------------------------------------------------------------

% Copyright (C) 2008-2009 
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
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

n_vett=[20; 40; 80; 100; 200; 300; 400; 500];  % PADUA POINTS DEGREES VECTOR.
xyrange=[-1,1,-1,1];                           % DEFINITION OF THE RECTANGLE.
number_test=10;                                % WE WILL MAKE "number_test" TO 
                                               % COMPUTE THE AVERAGE CPUTIME.
                                               
%-------------------------- SETTINGS END HERE ----------------------------------

for index=1:length(n_vett)
    
    n=n_vett(index);
    
    fprintf('\n \n \t [DEGREE]: %3.0f',n);
    
    % FFT METHOD.
    t(3)=cputime;
    for index=1:number_test
        wFFT = pdwtsFFT(n,xyrange);
    end
    t(4)=cputime;
    
    fprintf('\n \t [FFT][NUMBER OF TESTS]: %3.0f [AVERAGE CPUTIME]: %2.5f',...
        number_test,(t(4)-t(3))/number_test);
    
    % MM METHOD.
    t(1)=cputime;
    for index=1:number_test
        wMM = pdwtsMM(n,xyrange);
    end
    t(2)=cputime;
    
    fprintf('\n \t [MM] [NUMBER OF TESTS]: %3.0f [AVERAGE CPUTIME]: %2.5f',...
        number_test,(t(2)-t(1))/number_test);
    
end

fprintf('\n \n');