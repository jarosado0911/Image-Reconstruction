%This demonstrates the difference in adding up sorted vs nonsorted numbers
myContinue=0;
while myContinue~=1
    rge=input('Enter a range enclosed in brackets: ');
    x=input('How many random decimals? ');
    A=rge(1)+abs(rge(1)-rge(2))*rand(x,1);
    sumA=sum(A);
    sumSortA=sum(sort(A));
    diff=abs(sumA-sumSortA);
    sprintf('Using MatLab sum function: \n The sum of unsorted decimals = %d. \n The sum of the sorted decimes = %d. \n The difference in sums (sorted-unsorted) =%d.' ,sumA, sumSortA,diff)
    
    %Comp Sum of Unsorted numbers
    s1=compSum(A);
    diff1=abs(sumA-s1);
    diff2=abs(sumSortA-s1);
    sprintf('Using compensated the sum (unsorted) = %d \n The difference (MatLab Sum) from unsorted = %d and sorted = %d',s1, diff1,diff2)
    
    %Comp Sum of sorted numbers
    s2=compSum(sort(A));
    diff3=abs(sumA-s2);
    diff4=abs(sumSortA-s2);
    sprintf('Using compensated the sum (sorted) = %d \n The difference (MatLab Sum) from unsorted = %d and sorted = %d',s2, diff3,diff4)
    myContinue=input('Again? [0=yes/1=no]: ');
end