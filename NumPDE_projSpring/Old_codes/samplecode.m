for n=1:nT
    for i=2:npts-1
        for j=2:npts-1
            phi(i,j)=cfl*sqrt((phi(i,j)-phi(i,j-1))^2+(phi(i,j)-phi(i-1,j))^2)*g_DU(i,j)*kurv(i,j)...
                +(1+cfl*(Dxg_DU(i,j)+Dyg_DU(i,j)))*phi(i,j)-cfl*(Dxg_DU(i,j)*phi(i-1,j)+Dyg_DU(i,j)*phi(i,j-1));
            
            phi2(i,j)=(1+cfl*(Dxg_DU(i,j)+Dyg_DU(i,j)))*phi2(i,j)-cfl*(Dxg_DU(i,j)*phi2(i-1,j)+Dyg_DU(i,j)*phi2(i,j-1));
        end
    end