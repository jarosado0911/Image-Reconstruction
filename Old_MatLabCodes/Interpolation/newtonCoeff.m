%Newton Coefficients functions
function c_out = newtonCoeff(xlst,ylst)
flst=ylst;
listlenX=length(xlst);
        for k=1:listlenX
            temp1=flst(k);
                for n=k:listlenX-1
                    temp2=flst(n+1);
                    flst(n+1)=diffq(xlst(n-k+1),xlst(n+1),temp1,temp2);
                    temp1=temp2;
                end
        end
 c_out=flst;
end

