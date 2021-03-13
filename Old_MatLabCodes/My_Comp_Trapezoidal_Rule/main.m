%format long
clf
close all
clear
clc

myfun=@(x) (1+25*x.^2).^(-1);
mydfun=@(x) -(1+25*x.^2).^(-2).*(50.*x);
a=-1;
b=1;
nsamples=[1:10:110];
frm=ceil(length(nsamples)/2);

hfig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';


for j=1:length(nsamples)
    clf
    N=nsamples(j);
    h=abs(b-a)/N;
    xlst=linspace(a,b,N+1);
    flst=myfun(xlst);
    Iexact=integral(myfun,a,b);

    AreaSum=0;
    cumSum=[];
    subErrors=[];
    subTrueErr=[];
    subAreas=[];
for n=1:length(xlst)-1
    [tempA,tempE,tempTe]=TrapezoidArea(xlst(n),xlst(n+1),flst(n),flst(n+1),mydfun,myfun);
    subErrors=[subErrors,tempE];
    subAreas=[subAreas,tempA];
    subTrueErr=[subTrueErr,tempTe];
    AreaSum=AreaSum+tempA;
    cumSum=[cumSum,AreaSum];
end
    subErrors;
    AreaSum;
    compError=h^2/12*(mydfun(a)-mydfun(b));
    AreaSumComp=AreaSum+compError;
    diff=abs(Iexact-AreaSum);
    Ttab=table(xlst(1:N)',xlst(2:N+1)',subErrors', subTrueErr',subAreas',cumSum','VariableNames',{'xi','xf','subErrors','subTrueErr','subAreas','ACC_Area'});
    outA=[xlst(1:N)',xlst(2:N+1)',subErrors',subTrueErr',subAreas',cumSum'];
    fileID = fopen('Comp_Trap_Out.txt','a');

fprintf(fileID,'\n \n %s %s %s %s %s %s\r\n','xi','xf','subErrors','subTrueErr','subAreas','ACC_Area');
fprintf(fileID,'%d %d %d %d %d %d\r\n',outA);
fprintf(fileID,'The exact area of our function is = %d. \n The uncompensated trapezoid sum is = %d.\n The compensated trapezoid sum is = %d.\n The compensation error is = %d.\n The sum of the subErrors = %d.\n The sum of true SubErrors = %d.\n',Iexact,AreaSum,AreaSumComp,compError,sum(subErrors),sum(subTrueErr))

fprintf('The exact area of our function is = %d.\n The difference of exact and approx = %d using %u samples. \n The uncompensated trapezoid sum is = %d.\n The compensated trapezoid sum is = %d.\n The compensation error is = %d.\n The sum of the subErrors = %d.\n The sum of true SubErrors = %d.\n',Iexact,diff,nsamples(j),AreaSum,AreaSumComp,compError,sum(subErrors),sum(subTrueErr))

%subplot(2,frm,j)
hold on
fplot(myfun,[a,b],'b')
stem(xlst,flst,'r')

for k=1:1:length(xlst)-1
    rx = [xlst(k) xlst(k) xlst(k+1) xlst(k+1)];
    ry = [flst(k) 0 flst(k+1) 0];
    k = convhull(rx, ry);
    fill (rx(k), ry(k), 'g','facealpha', 0.23); 
end

drawnow 
     % Capture the plot as an image 
      frame = getframe(hfig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
fclose(fileID);