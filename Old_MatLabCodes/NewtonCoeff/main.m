prompt1='Please enter in a list of x-values using brackets! ';
prompt2='Please enter in a list of corresponding y-values using brackets!: ';
prompt3='Resolution, enter in a posiive integer? ';
promptend='Again? (0=yes/1=no) ';
again=0;
while again ~=1
    xlst=input(prompt1);
    flst=input(prompt2);
    res=input(prompt3);
    maxX=max(xlst);
    minX=min(xlst);
    myX=linspace(minX,maxX,res);
    listlenX= length(xlst);
    if length(flst)== listlenX
        disp('Your lists are okay!');
        disp('Here are your coefficients!')
        myCo=newtonCoeff(xlst,flst);
        disp(myCo);
        myY=newtonP(myX,xlst,myCo);
        scatter(myX,myY,20);
        hold on
        scatter(xlst,flst,140,'d','filled','MarkerFaceColor',[0 .7 .7]);
    else
        disp('Your lists need to be same length!');
    end
    again=input(promptend);
    hold off
end

%{
xlst=[1,2,3,4];
flst=[5,8,5,1];

flst
for k=1:4
temp1=flst(k);
for n=k:3
    temp2=flst(n+1);
    flst(n+1)=diffq(xlst(n-k+1),xlst(n+1),temp1,temp2);
    temp1=temp2;
end
end
flst
%}
