  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHORS:

       Marco Caliari
       University of Verona, Italy
       E-mail: marco.caliari@univr.it

       Stefano de Marchi, Alvise Sommariva, Marco Vianello
       University of Padua, Italy
       E-mail: demarchi@math.unipd.it, alvise@math.unipd.it,
               marcov@math.unipd.it

   REFERENCE:

    -  Padua2DM: fast interpolation and cubature at the Padua points in
       Matlab/Octave
       NUMERICAL ALGORITHMS, 56 (2011), PP. 45-60

   SOFTWARE REVISION DATE:

       V1.0, January 2010

   SOFTWARE LANGUAGE:

       MATLAB 7.6 and Octave 3.2


================================================================================
SOFTWARE
================================================================================
This software provides all the functions needed to perform interpolation and
cubature at the Padua points, together with the functions and the demos used
in the paper.

================================================================================
PACKAGE
================================================================================

The directory contains the following files

README.txt              : this file
pdint.m                 : main function for interpolation at the Padua points
pdcub.m                 : main function for cubature at the Padua points
pdpts.m                 : function for the computation of the Padua points
pdcfsFFT.m              : function for the computation of the interpolation
                          coefficients by FFT (recommended)
pdcfsMM.m               : function for the computation of the interpolation
                          coefficients by matrix multiplications
pdval.m                 : function for the evaluation of the interpolation
                          polynomial
pdwtsFFT.m              : function for the computation of the cubature
                          weights by FFT
pdwtsMM.m               : function for the computation of the cubature
                          weights by matrix multiplications (recommended)
funct.m                 : function containing some test functions
demo_pdint.m            : demo script for pdint
demo_cputime_pdint.m    : demo script for the computation of CPU time for
                          interpolation
demo_errors_pdint.m     : demo script for the comparison of interpolation with
                          coefficients computed by FFT or by matrix
                          multiplications
demo_pdcub              : demo script for pdcub
demo_cputime_pdcub.m    : demo script for the computation of CPU time for
                          cubature
demo_errors_pdcub.m     : demo script for the comparison of cubature with
                          weights computed by FFT or by matrix multiplications
demo_errors_pdcub_gl.m  : demo script for the comparison of different cubature
                          formulas
cubature_square.m       : function for the computation of some cubature
                          formulas for the square
omelyan_solovyan_rule.m : function for the computation of Omelyan-Solovyan
                          cubature points and weights
Contents.m              : Contents file for Matlab

==============
HOW TO INSTALL
==============

Unpack Padua2DM.tgz. The directory Padua2DM will be created with all
the needed files inside. Add the directory Padua2DM to the path in
Matlab/Octave or change the current working directory to Padua2DM.

================================================================================
TEST
================================================================================

The core files pdint.m, pdcub.m, pdpts.m, pdcfsFFT.m, pdcfsMM.m, pdval.m,
pdwtsFFT.m, pdwtsMM.m, testfunct.m contain the Octave test facility. It is
possibile to enter, for instance,

octave:1> test pdint

to check if the function performs as expected and passes all the built-in tests.
