function J=jGaussian(x,y,var)

[X,Y]=meshgrid(x,y);
J=zeros(length(x),length(y));

for i=1:length(x)
    for j=1:length(y)
        J(i,j)=(1/(4*pi*var))*exp((X(i,j)^2+Y(i,j)^2)/(4*var));
    end
end
end

