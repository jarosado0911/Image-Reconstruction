function y_out = newtonP(xin,xlst,nCoeff)
sum=0;
for k=1:length(xlst)
    prod=nCoeff(k);
    for j=1:k-1
        prod=prod.*(xin-xlst(j));
    end
    sum=sum+prod;
end
y_out=sum;
end

