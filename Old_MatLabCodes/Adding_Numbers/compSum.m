function s = compSum(lst)
    s=0;
    c=0;
    for j=1:length(lst)
        temp=s;
        y=lst(j)+c;
        s=temp+y;
        c=(temp-s)+y;
    end
end

