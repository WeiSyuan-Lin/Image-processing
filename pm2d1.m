% 2D Perona Malik solver 
% physcial parameters

clc
clear all;
% numerical parameters
nu=1;
N=40;
M=40;


[N,M]=size(I2);

dx=1/N;
dy=1/M;
dt=0.1*dx^2/nu;
lambda1=dt*nu/(dx)^2;
lambda2=dt*nu/(dy)^2;
maxiter=1/dt;
x=dx*(1/2:1:N-1/2);
y=dy*(1/2:1:M-1/2);
[X, Y]=meshgrid(x,y);

% initial data

U=zeros(N,M);
U=cos(pi.*X).*cos(pi.*Y);

U=U';

C=zeros(N+2,M+2);
K=3;
%coef=@(x) nu;
coef=@(x) 1./sqrt(1+x);
%coef=@(x) 1./(1+x./K);
%coef=@(x)  exp(-x./K);
%

% forcing term
force=@(xt,yt,t) -exp(-t).*cos(pi.*xt).*cos(pi.*yt) ...
                +2*pi.^2.*exp(-t).*cos(pi.*xt).*cos(pi.*yt) ...
                ./sqrt(1+pi.^2.*exp(-2.*t).*(sin(pi.*xt).^2.*cos(pi.*yt).^2 +cos(pi.*xt).^2.*sin(pi.*yt).^2)) ...
                -pi.^4.*exp(-3.*t).*(cos(2*pi.*xt).*sin(pi.*yt).^2+sin(pi.*xt).^2.*cos(2*pi.*yt)).*cos(pi.*xt).*cos(pi.*yt) ...
               ./(1+pi.^2.*exp(-2.*t).*(sin(pi.*xt).^2.*cos(pi.*yt).^2 +cos(pi.*xt).^2.*sin(pi.*yt).^2)).^(3/2);
            
           
            

 for k=1:maxiter
 U_old=U;
 
 %update coefficient
 % update i:3:N, j:3:M
 Dx=(U(3:N,2:M-1)-U(1:N-2,2:M-1))./(2*dx);
 Dy=(U(2:N-1,3:M)-U(2:N-1,1:M-2))./(2*dy);
 C(3:N,3:M)=coef(Dx.^2+Dy.^2);
 
 % update i:2, j:3:M
 Dx=(U(2,2:M-1)-U(1,2:M-1))./(2*dx);
 Dy=(U(1,3:M)-U(1,1:M-2))./(2*dy);
 C(2,3:M)=coef(Dx.^2+Dy.^2);
 
 % update i:N+1, j:3:M
 Dx=(U(N,2:M-1)-U(N-1,2:M-1))./(2*dx);
 Dy=(U(N,3:M)-U(N,1:M-2))./(2*dy);
 C(N+1,3:M)=coef(Dx.^2+Dy.^2);
 
 % update i:3:N, j:2
 Dx=(U(3:N,1)-U(1:N-2,1))./(2*dx);
 Dy=(U(2:N-1,2)-U(2:N-1,1))./(2*dy);
 C(3:N,2)=coef(Dx.^2+Dy.^2);
 
 % update i:3:N, j:M+1
 Dx=(U(3:N,M)-U(1:N-2,M))./(2*dx);
 Dy=(U(2:N-1,M)-U(2:N-1,M-1))./(2*dy);
 C(3:N,M+1)=coef(Dx.^2+Dy.^2);
 
 % update i:2, j:2
 Dx=(U(2,1)-U(1,1))./(2*dx);
 Dy=(U(1,2)-U(1,1))./(2*dy);
 C(2,2)=coef(Dx.^2+Dy.^2);
 
 % update i:2, j:M+1
 Dx=(U(2,M)-U(1,M))./(2*dx);
 Dy=(U(1,M)-U(1,M-1))./(2*dy);
 C(2,M+1)=coef(Dx.^2+Dy.^2);
 
 % update i:N+1, j:M+1
 Dx=(U(N,M)-U(N-1,M))./(2*dx);
 Dy=(U(N,M)-U(N,M-1))./(2*dy);
 C(N+1,M+1)=coef(Dx.^2+Dy.^2);
 
 % update i:N+1, j:2
 Dx=(U(N,1)-U(N-1,1))./(2*dx);
 Dy=(U(N,2)-U(N,1))./(2*dy);
 C(N+1,2)=coef(Dx.^2+Dy.^2);
 
  % update i:1, j:3:M
 Dx=(U(1,2:M-1)-U(2,2:M-1))./(2*dx);
 Dy=(U(1,3:M)-U(1,1:M-2))./(2*dy);
 C(1,3:M)=coef(Dx.^2+Dy.^2);
 
  % update i:N+2, j:3:M
 Dx=(U(N-1,2:M-1)-U(N,2:M-1))./(2*dx);
 Dy=(U(N,3:M)-U(N,1:M-2))./(2*dy);
 C(N+2,3:M)=coef(Dx.^2+Dy.^2);
 
 % update i:3:N, j:1
 Dx=(U(3:N,1)-U(1:N-2,1))./(2*dx);
 Dy=(U(2:N-1,1)-U(2:N-1,2))./(2*dy);
 C(3:N,1)=coef(Dx.^2+Dy.^2);
 
 % update i:3:N, j:M+2
 Dx=(U(3:N,M)-U(1:N-2,M))./(2*dx);
 Dy=(U(2:N-1,M-1)-U(2:N-1,M))./(2*dy);
 C(3:N,M+2)=coef(Dx.^2+Dy.^2);
 
 U_old(2:N-1, 2:M-1)=U(2:N-1, 2:M-1) ...
  +0.5*lambda1*((C(4:N+1,3:M)+C(3:N,3:M)).*U(3:N,2:M-1)...
  -(2*C(3:N,3:M)+C(4:N+1,3:M)+C(2:N-1,3:M)).*U(2:N-1,2:M-1)...
              +(C(2:N-1,3:M)+C(3:N,3:M)).*U(1:N-2,2:M-1)) ...
  +0.5*lambda2*((C(3:N,4:M+1)+C(3:N,3:M)).*U(2:N-1,3:M)...
  -(C(3:N,4:M+1)+2*C(3:N,3:M)+C(3:N,2:M-1)).*U(2:N-1,2:M-1)...
    +(C(3:N,2:M-1)+C(3:N,3:M)).*U(2:N-1,1:M-2));

U_old(2:N-1, 2:M-1)=U_old(2:N-1, 2:M-1) + dt*force(X(2:N-1,2:M-1), Y(2:N-1,2:M-1), (k-1)*dt);

% i=1 update
U_old(1, 2:M-1)=U(1, 2:M-1) ...
  +0.5*lambda1*((C(3,3:M)+C(2,3:M)).*U(2,2:M-1)...
  -(2*C(2,3:M)+C(3,3:M)+C(1,3:M)).*U(1,2:M-1)...
              +(C(1,3:M)+C(2,3:M)).*U(1,2:M-1)) ...
  +0.5*lambda2*((C(2,4:M+1)+C(2,3:M)).*U(1,3:M)...
  -(C(2,4:M+1)+2*C(2,3:M)+C(2,2:M-1)).*U(1,2:M-1)...
    +(C(2,2:M-1)+C(2,3:M)).*U(1,1:M-2));

U_old(1, 2:M-1)=U_old(1, 2:M-1) + dt*force(X(1,2:M-1), Y(1,2:M-1), (k-1)*dt);

% i=N update
 U_old(N, 2:M-1)=U(N, 2:M-1) ...
  +0.5*lambda1*((C(N+2,3:M)+C(N+1,3:M)).*U(N,2:M-1)...
  -(2*C(N+1,3:M)+C(N+2,3:M)+C(N,3:M)).*U(N,2:M-1)...
              +(C(N,3:M)+C(N+1,3:M)).*U(N-1,2:M-1)) ...
  +0.5*lambda2*((C(N+1,4:M+1)+C(N+1,3:M)).*U(N,3:M)...
  -(C(N+1,4:M+1)+2*C(N+1,3:M)+C(N+1,2:M-1)).*U(N,2:M-1)...
    +(C(N+1,2:M-1)+C(N+1,3:M)).*U(N,1:M-2));

U_old(N, 2:M-1)=U_old(N, 2:M-1) + dt*force(X(N,2:M-1), Y(N,2:M-1), (k-1)*dt);
 
% j=1 update
 U_old(2:N-1, 1)=U(2:N-1, 1) ...
  +0.5*lambda1*((C(4:N+1,2)+C(3:N,2)).*U(3:N,1)...
  -(2*C(3:N,2)+C(4:N+1,2)+C(2:N-1,2)).*U(2:N-1,1)...
              +(C(2:N-1,2)+C(3:N,2)).*U(1:N-2,1)) ...
  +0.5*lambda2*((C(3:N,3)+C(3:N,2)).*U(2:N-1,2)...
  -(C(3:N,3)+2*C(3:N,2)+C(3:N,1)).*U(2:N-1,1)...
    +(C(3:N,1)+C(3:N,2)).*U(2:N-1,1));

 U_old(2:N-1, 1)=U_old(2:N-1, 1) + dt*force(X(2:N-1,1), Y(2:N-1,1), (k-1)*dt);


% j=M update
 U_old(2:N-1, M)=U(2:N-1, M) ...
  +0.5*lambda1*((C(4:N+1,M+1)+C(3:N,M+1)).*U(3:N,M)...
  -(2*C(3:N,M+1)+C(4:N+1,M+1)+C(2:N-1,M+1)).*U(2:N-1,M)...
              +(C(2:N-1,M+1)+C(3:N,M+1)).*U(1:N-2,M)) ...
  +0.5*lambda2*((C(3:N,M+2)+C(3:N,M+1)).*U(2:N-1,M)...
  -(C(3:N,M+2)+2*C(3:N,M+1)+C(3:N,M)).*U(2:N-1,M)...
    +(C(3:N,M)+C(3:N,M+1)).*U(2:N-1,M-1));

U_old(2:N-1, M)=U_old(2:N-1, M) + dt*force(X(2:N-1,M), Y(2:N-1,M), (k-1)*dt);


 % i=1, j=1 update
U_old(1, 1)=U(1, 1) ...
  +0.5*lambda1*((C(3,2)+C(2,2)).*U(2,1)...
  -(2*C(2,2)+C(3,2)+C(1,2)).*U(1,1)...
              +(C(1,2)+C(2,2)).*U(1,1)) ...
  +0.5*lambda2*((C(2,3)+C(2,2)).*U(1,2)...
  -(C(2,3)+2*C(2,2)+C(2,1)).*U(1,1)...
    +(C(2,2)+C(2,1)).*U(1,1));

U_old(1, 1)=U_old(1, 1) + dt*force(X(1,1), Y(1,1), (k-1)*dt);

% i=N, j=1 update
 U_old(N, 1)=U(N, 1) ...
  +0.5*lambda1*((C(N+2,2)+C(N+1,2)).*U(N,1)...
  -(2*C(N+1,2)+C(N+2,2)+C(N,2)).*U(N,1)...
              +(C(N,2)+C(N+1,2)).*U(N-1,1)) ...
  +0.5*lambda2*((C(N+1,3)+C(N+1,2)).*U(N,2)...
  -(C(N+1,3)+2*C(N+1,2)+C(N+1,1)).*U(N,1)...
    +(C(N+1,2)+C(N+1,1)).*U(N,1));
 
U_old(N, 1)=U_old(N, 1) + dt*force(X(N,1), Y(N,1), (k-1)*dt);  

% j=M, i=1 update
 U_old(1, M)=U(1, M) ...
  +0.5*lambda1*((C(3,M+1)+C(2,M+1)).*U(2,M)...
  -(2*C(2,M+1)+C(3,M+1)+C(1,M+1)).*U(1,M)...
              +(C(1,M+1)+C(2,M+1)).*U(1,M)) ...
  +0.5*lambda2*((C(2,M+2)+C(2,M+1)).*U(1,M)...
  -(C(2,M+2)+2*C(2,M+1)+C(2,M)).*U(1,M)...
    +(C(2,M+1)+C(2,M)).*U(1,M-1));
U_old(1, M)=U_old(1, M) + dt*force(X(1,M), Y(1,M), (k-1)*dt);

% j=M, i=N update
 U_old(N, M)=U(N, M) ...
  +0.5*lambda1*((C(N+2,M+1)+C(N+1,M+1)).*U(N,M)...
  -(2*C(N+1,M+1)+C(N+2,M+1)+C(N,M+1)).*U(N,M)...
              +(C(N,M+1)+C(N+1,M+1)).*U(N-1,M)) ...
  +0.5*lambda2*((C(N+1,M+2)+C(N+1,M+1)).*U(N,M)...
  -(C(N+1,M+2)+2*C(N+1,M+1)+C(N+1,M)).*U(N,M)...
    +(C(N+1,M+1)+C(N+1,M)).*U(N,M-1));
U_old(N, M)=U_old(N, M) + dt*force(X(N,M), Y(N,M), (k-1)*dt);

U=U_old;
 end
 
 %error 
 
 % exact solution
 
 U_exact=zeros(N,M);
 U_exact=exp(-dt*maxiter)*cos(pi.*X).*cos(pi.*Y);
 U_exact=U_exact';
 
 norm(U_exact-U)*sqrt(dx*dy)
 surf(U)%
                              
                              


% 