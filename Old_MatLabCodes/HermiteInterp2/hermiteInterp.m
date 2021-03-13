%HERMITE INTERPOLATION FUNCTION
function coeffs = hermiteInterp(Xlst,Ylst,YDlst,YDDlst,polDeg)
Nmat=zeros(polDeg+1);
Nmat(:,1)=Xlst';
Nmat(:,2)=Ylst';
Nmat(:,3)=YDlst';
Nmat(:,4)=YDDlst';
Nmat;
for j=1:polDeg
    for k=1:length(Xlst)-j
        if Nmat(k,1) ~= Nmat(k+j,1)
            Nmat(k,j+2)=(Nmat(k,j+1)-Nmat(k+1,j+1))/((Nmat(k,1)-Nmat(k+j,1)));
        end
    end
end
coeffs=Nmat(1,2:polDeg+2);
end

