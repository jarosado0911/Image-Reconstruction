function y_out = polyNum(xin,coeffs)
sum=0;
for k=1:length(coeffs)
    prod=coeffs(k);
    for j=1:k-1
        prod=prod.*xin;
    end
    sum=sum+prod;
end
y_out=sum;
end

