function U=gauss2d(raster,x,nx,sigma)

u = raster;

% Problem parameters
tf = sigma; % final time
nt = 10; % number of time steps
% Iniatialization
k = tf/nt; % time step
h = x(2)-x(1); % desired mesh resolution

U = reshape(u,[],1); % initial conditions

% Matrix generation
K1D = spdiags(ones(nx,1)*[1 -2 1],-1:1,nx,nx)/h^2; % 1d Poisson matrix
I1D = speye(size(K1D)); % 1d identity matrix
K2D = kron(K1D,I1D)+kron(I1D,K1D); % 2d Poisson matrix
I2D = speye(size(K2D)); % 2d identity matrix
A = I2D-k/2*K2D; % LHS matrix for Crank-Nicolson
B = I2D+k/2*K2D; % RHS matrix for Crank-Nicolson

ndpts = floor(3);
s = mit18086_stencil_stability((0:1:ndpts),1);
% Time loop
for j = 0:nt
    if j>0 % update step
        tic
        U = A\(B*U);
    end
    
    % for neumann b.c.
    U = weightedsum(U,s,nx);
end
U = reshape(U,nx,nx);
end

function U=weightedsum(U,s,nx)
U = reshape(U,nx,nx);
w=s(2:end);
U(1,:)=0*U(1,:);
U(:,1)=0*U(:,1);
U(end,:)=U(1,:);
U(:,end)=U(:,1);

for j=1:length(w)
 U(1,:)=U(1,:)-(w(j)/s(1))*U(j+1,:);
 U(:,1)=U(:,1)-(w(j)/s(1))*U(:,j+1);
 U(end,:)=U(end,:)-(w(j)/s(1))*U(nx-j,:);
 U(:,end)=U(:,end)-(w(j)/s(1))*U(:,nx-j);
end

U = reshape(U,[],1);
end