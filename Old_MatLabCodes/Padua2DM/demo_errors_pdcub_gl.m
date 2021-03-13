%-------------------------------------------------------------------------------
% This demo computes the cubature of a function "f" on a rectangle by an 
% algebraic formula of degree "n" on Padua points as well as by other formulas
% (tensorial Gauss-Legendre-Lobatto, tensorial Gauss-Legendre, tensorial 
% Clenshaw-Curtis, Osmolyan-Solovyan minimal rule) on a rectangle. 
% A comparison of relative errors between FFT method and MM method is made.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdpts
% 2. pdwtsFFT
% 3. cubature_square
% 4. omelyan_solovyan_rule
% 5. funct
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

n_vett=2:1:30;                            % PADUA POINTS DEGREES VECTOR.
xyrange=[-1,1,-1,1];                      % DEFINITION OF THE RECTANGLE.
                                          
n_exact=200;                              % COMPUTING APPROX. OF EXACT INTEGRAL.

%-------------------------- SETTINGS END HERE ----------------------------------

PadL_exact = pdpts(n_exact,xyrange);
fPadL_exact=feval(@testfunct,PadL_exact(:,1),PadL_exact(:,2));
w_exact = pdwtsFFT(n_exact,xyrange);
I_exact = w_exact'*fPadL_exact;                                           
        
cubOS=[]; number_ptsOS=[]; relerrOS=[];   % STORAGE OF OMELYAN-SOLOVYAN RULE.

for index=1:length(n_vett)
    
    n=n_vett(index);
    
    fprintf('\n \n \t [DEGREE]: %3.0f',n);
    PadL = pdpts(n,xyrange);
    fPadL=feval(@testfunct,PadL(:,1),PadL(:,2));
    
    % PADUA-FFT METHOD.
    wFFT = pdwtsFFT(n,xyrange);
    cubFFT(index) = wFFT'*fPadL;
    number_ptsFFT(index)=size(wFFT,1)*size(wFFT,2);
    Iden=(1-abs(sign(I_exact)))+I_exact;
    relerrFFT(index)=abs((cubFFT(index)-I_exact)/Iden);
    fprintf('\n \n \t [DEG]: %4.0f [PTS.PD.]: %5.0f [REL.ERR.PD.]: %2.5e', ...
    n,number_ptsFFT(index),relerrFFT(index));

    % TENSORIAL CLENSHAW-CURTIS METHOD.
    [xv,yv,wv]=cubature_square(n,xyrange(1),xyrange(2),xyrange(3),xyrange(4),3);
    fv=feval(@testfunct,xv,yv);
    cubCC(index)=wv'*fv;
    number_ptsCC(index)=size(wv,1)*size(wv,2);
    relerrCC(index)=abs((cubCC(index)-I_exact)/Iden);
    fprintf('\n \t [DEG]: %4.0f [PTS.CC.]: %5.0f [REL.ERR.CC.]: %2.5e', ...
    n,number_ptsCC(index),relerrCC(index));

    % TENSORIAL GAUSS-LEGENDRE METHOD.
    [xv,yv,wv]=cubature_square(n,xyrange(1),xyrange(2),xyrange(3),xyrange(4),4);
    fv=feval(@testfunct,xv,yv);
    cubGL(index)=wv'*fv;
    number_ptsGL(index)=size(wv,1)*size(wv,2);
    relerrGL(index)=abs((cubGL(index)-I_exact)/Iden);
    fprintf('\n \t [DEG]: %4.0f [PTS.GL.]: %5.0f [REL.ERR.GL.]: %2.5e', ...
    n,number_ptsGL(index),relerrGL(index));

    % TENSORIAL GAUSS-LEGENDRE-LOBATTO METHOD.
    [xv,yv,wv]=cubature_square(n,xyrange(1),xyrange(2),xyrange(3),xyrange(4),5);
    fv=feval(@testfunct,xv,yv);
    cubGLL(index)=wv'*fv;
    number_ptsGLL(index)=size(wv,1)*size(wv,2);
    relerrGLL(index)=abs((cubGLL(index)-I_exact)/Iden);
    fprintf('\n \t [DEG]: %4.0f [PTS.GLL.]: %5.0f [REL.ERR.GLL.]: %2.5e', ...
    n,number_ptsGLL(index),relerrGLL(index));

    % OMELYAN-SOLOVYAN RULE.
    if (n == 15) | (n == 17) | (n == 19) | (n == 21) | (n == 23)
       [xv,yv,wv]=omelyan_solovyan_rule(n);
       fv=feval(@testfunct,xv,yv);
       cubOS_loc=wv'*fv; 
       cubOS=[cubOS; cubOS_loc]; 
       
       number_ptsOS_loc=size(wv,1)*size(wv,2); 
       number_ptsOS=[number_ptsOS; number_ptsOS_loc];
       
       relerrOS_loc=abs((cubOS_loc-I_exact)/Iden);
       relerrOS=[relerrOS; relerrOS_loc];
       fprintf('\n \t [DEG]: %4.0f [PTS.OS.]: %5.0f [REL.ERR.OS.]: %2.5e', ...
                 n,number_ptsOS_loc,relerrOS_loc);
   end
   
end

% PLOT (PTS - ERRORS OF PADUA FFT) vs. 
% (PTS - ERRORS OF TENSORIAL GAUSS-LEGENDRE)

[disp_ptsFFT, indexFFT]=find(number_ptsFFT <= 500);
[disp_ptsCC, indexCC]=find(number_ptsCC <= 500);
[disp_ptsGL, indexGL]=find(number_ptsGL <= 500);
[disp_ptsGLL, indexGLL]=find(number_ptsGLL <= 500);
[disp_ptsOS, indexOS]=find(number_ptsOS <= 500);

if length(indexOS > 0)
    semilogy(number_ptsFFT(indexFFT),relerrFFT(indexFFT),'g-d', ...
        number_ptsCC(indexCC),relerrCC(indexCC),'b-x',number_ptsGL(indexGL), ...
        relerrGL(indexGL),'m-o',number_ptsGLL(indexGLL),relerrGLL(indexGLL), ...
        'r-*',number_ptsOS(disp_ptsOS),relerrOS(disp_ptsOS),'c-s');
    legend('Non tens. CC Padua pts.','Tens. CC','Tens. GL','Tens. GLL', ...
        'Non tens. OS pts.');
else
    semilogy(number_ptsFFT(indexFFT),relerrFFT(indexFFT),'g-d', ...
        number_ptsCC(indexCC),relerrCC(indexCC),'b-x',number_ptsGL(indexGL), ...
        relerrGL(indexGL),'m-o',number_ptsGLL(indexGLL),relerrGLL(indexGLL), ...
        'r-*');
    legend('Non tens. CC Padua pts.','Tens. CC','Tens. GL','Tens. GLL');
end
xlabel('Number of points');
ylabel('Relative cubature error');

fprintf('\n \n');