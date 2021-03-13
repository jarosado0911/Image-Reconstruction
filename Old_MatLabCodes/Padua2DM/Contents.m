% Padua2DM: fast interpolation and cubature at the Padua points 
% in Matlab/Octave. 
%
% This software provides all the functions needed to perform
% interpolation and cubature at the Padua points, together with
% the functions and the demos used in the paper.
%
% pdint     - The main function for interpolation at the Padua points.
% pdcub     - The main function for cubature at the Padua points.
% pdpts     - Function for the computation of the Padua points
% pdcfsFFT  - Function for the computation of the interpolation coefficients 
%             by FFT (recommended).
% pdcfsMM   - Function for the computation of the interpolation coefficients 
%             by matrix multiplications
% pdval     - Function for the evaluation of the interpolation polynomial
% pdwtsFFT  - Function for the computation of the cubature weights by FFT
% pdwtsMM   - Function for the computation of the cubature weights 
%             by matrix multiplications (recommended)
% testfunct - Function containing some test functions
%
% The package provides also several demo scripts.
% For more details concerning input and output parameters and the usage of 
% single functions, see the corresponding help.
%
% Authors:
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
