
N=6;
M=8;
dx=1/N;
dy=1/M;
mu=1;
dt=1;
lamb1=(dt*mu)/(dx)^2; lamb2=(dt*mu)/(dy)^2;

U=zeros(N,M);
U=cos(pi.*X).*cos(pi.*Y);

d1=[1-lam1-lam2 repmat(1-2*lam1-lam2,[1,N-2]) 1-lam1-lam2];
d2=[1-lam1-2*lam2 repmat(1-2*lam1-2*lam2,[1,N-2]) 1-lam1-2*lam2];

%A=diag([d1 repmat(d2,[1,N-2]) d1])+lamb1*diag(ones(N*M-1,1),1)...
 %   +lamb1*diag(ones(N*M-1,1),-1)+lamb2*diag(ones(N*M-N,1),N)...
  %  +lamb2*diag(ones(N*M-N,1),-N);
  
  
 for k=1:10
 
 U_old=U;
 U_old(2:N-1,2:M-1)=U(2:N-1,2:M-1)...
                        +lamb1*U(3:N,2:M-1)-2*U(2:N-1,2:M-1)+U(1:N-2,2:M-1)...
                        +lamb2*U(2:N-1,3:M)-2*U(2:N-1,2:M-1)+U(2:N-1,1:M-2);

 U_old(1,2:M-1)=U(1,2:M-1)...
                        +lamb1*U(2:N,2:M-1)-2*U(1:N-1,2:M-1)+U(N-1,2:M-1)...
                        +lamb2*U(N,3:M)-2*U(N,2:M-1)+U(N,1:M-2);
        
U_old(N,2:M-1)=U(N,2:M-1)...
                        +lamb1*U(2:N,2:M-1)-2*U(1:N-1,2:M-1)+U(N-1,2:M-1)...
                        +lamb2*U(N,3:M)-2*U(N,2:M-1)+U(N,1:M-2);

                    
 U_old(2:N-1,1)=U(2:N-1,1)...
                        +lamb1*U(3:N,2:M-1)-2*U(2:N-1,2:M-1)+U(1:N-2,2:M-1)...
                        +lamb2*U(2:N-1,3:M)-2*U(2:N-1,2:M-1)+U(2:N-1,1:M-2);

                    
 
 end

