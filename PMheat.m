% 2D heat solver 

%  clc
%  clear all;
% % numerical parameters
% physcial parameters

I=imread('f01.jpg');
I1=rgb2gray(I);
I2=double(I1);

J1=imnoise(I1,'salt & pepper',0.02);

J2=double(J1);


nu=1;
[N,M]=size(I2);

dx=1/N;
dy=1/M;
dt=0.1*dx^2/nu;
lambda1=dt*nu/(dx)^2;
lambda2=dt*nu/(dy)^2;
maxiter=800;%1/dt;
x=dx*(1/2:1:N-1/2);
y=dy*(1/2:1:M-1/2);
[X, Y]=meshgrid(x,y);

% initial data

% U=zeros(N,M);
% U=cos(pi.*X).*cos(pi.*Y);
% U=U';

U=zeros(N,M);
U=J2;

C=zeros(N+2,M+2);
K=10;
%coef=@(x)1./sqrt(1+x);
coef=@(x)1./sqrt(1+x./K);
%coef=@(x)exp(-x./K); 

H=fspecial('gaussian',10,50);


 for k=1:maxiter
 U_old=U;
 Us=conv2(U,H,'same');
 
 % update coef i=3:N j=3:M
 Dx=(Us(3:N,2:M-1)-Us(1:N-2,2:M-1))./(2*dx);
 Dy=(Us(2:N-1,3:M)-Us(2:N-1,1:M-2))./(2*dy);
 C(3:N,3:M)=coef(Dx.^2+Dy.^2);

 %i=2 j=3:M
 Dx=(Us(2,2:M-1)-Us(1,2:M-1))./(2*dx);
 Dy=(Us(1,3:M)-Us(1,1:M-2))./(2*dy);
 C(2,3:M)=coef(Dx.^2+Dy.^2);
 
 %i=N+1 j=3:M
 Dx=(Us(N,2:M-1)-Us(N-1,2:M-1))./(2*dx);
 Dy=(Us(N,3:M)-Us(N,1:M-2))./(2*dy);
 C(N+1,3:M)=coef(Dx.^2+Dy.^2);
 
 %j=2 i=3:N
 Dx=(Us(3:N,1)-Us(1:N-2,1))./(2*dx);
 Dy=(Us(2:N-1,2)-Us(2:N-1,1))./(2*dy);
 C(3:N,2)=coef(Dx.^2+Dy.^2);
 
 %j=M+1 i=3:N
 Dx=(Us(3:N,M)-Us(1:N-2,M))./(2*dx);
 Dy=(Us(2:N-1,M)-Us(2:N-1,M-1))./(2*dy);
 C(3:N,M+1)=coef(Dx.^2+Dy.^2);
 
 %i=2 j=2
 Dx=(Us(2,1)-Us(1,1))./(2*dx);
 Dy=(Us(1,2)-Us(1,1))./(2*dy);
 C(2,2)=coef(Dx.^2+Dy.^2);
 
 %i=2 j=M+1
 Dx=(Us(2,M)-Us(1,M))./(2*dx);
 Dy=(Us(1,M)-Us(1,M-1))./(2*dy);
 C(2,M+1)=coef(Dx.^2+Dy.^2);
 

 %i=N+1 j=M+1
 Dx=(Us(N,M)-Us(N-1,M))./(2*dx);
 Dy=(Us(N,M)-Us(N,M-1))./(2*dy);
 C(N+1,M+1)=coef(Dx.^2+Dy.^2);

 %i=N+1 j=2
 Dx=(Us(N,1)-Us(N-1,1))./(2*dx);
 Dy=(Us(N,2)-Us(N,1))./(2*dy);
 C(N+1,2)=coef(Dx.^2+Dy.^2);

 %i=1 j=3:M
 Dx=(Us(1,2:M-1)-Us(2,2:M-1))./(2*dx);
 Dy=(Us(1,3:M)-Us(1,1:M-2))./(2*dy);
 C(1,3:M)=coef(Dx.^2+Dy.^2);
 
 %i=N+2 j=3:M
 Dx=(Us(N-1,2:M-1)-Us(N,2:M-1))./(2*dx);
 Dy=(Us(N,3:M)-Us(N,1:M-2))./(2*dy);
 C(N+2,3:M)=coef(Dx.^2+Dy.^2);
 
 %j=1 i=3:N
 Dx=(Us(3:N,1)-Us(1:N-2,1))./(2*dx);
 Dy=(Us(2:N-1,1)-Us(2:N-1,2))./(2*dy);
 C(3:N,1)=coef(Dx.^2+Dy.^2);
 
 %j=M+2 i=3:N
 Dx=(Us(3:N,M)-Us(1:N-2,M))./(2*dx);
 Dy=(Us(2:N-1,M-1)-Us(2:N-1,M))./(2*dy);
 C(3:N,M+2)=coef(Dx.^2+Dy.^2);

  U_old(2:N-1,2:M-1)=U(2:N-1,2:M-1)...
      +0.5*lambda1*((C(4:N+1,3:M)+C(3:N,3:M)).*U(3:N,2:M-1)...
      -(C(4:N+1,3:M)+2*C(3:N,3:M)+C(2:N-1,3:M)).*U(2:N-1,2:M-1)...
      +(C(2:N-1,3:M)+C(3:N,3:M)).*U(1:N-2,2:M-1))...
     +0.5*lambda2*((C(3:N,4:M+1)+C(3:N,3:M)).*U(2:N-1,3:M)...
  -(C(3:N,4:M+1)+2*C(3:N,3:M)+C(3:N,2:M-1)).*U(2:N-1,2:M-1)...
  +(C(3:N,2:M-1)+C(3:N,3:M)).*U(2:N-1,1:M-2));

% i=1 update
  U_old(1, 2:M-1)=U(1, 2:M-1) ...
  +0.5*lambda1*((C(3,3:M)+C(2,3:M)).*U(2,2:M-1)...
  -(C(3,3:M)+2*C(2,3:M)+C(1,3:M)).*U(1,2:M-1)...
  +(C(2,3:M)+C(1,3:M)).*U(1,2:M-1)) ...
  +0.5*lambda2*((C(2,4:M+1)+C(2,3:M)).*U(1,3:M)...
  -(C(3,4:M+1)+2*C(2,3:M)+C(2,2:M-1)).*U(1,2:M-1)...
  +(C(2,2:M-1)+C(2,3:M)).*U(1,1:M-2));

% i=N update
  U_old(N, 2:M-1)=U(N, 2:M-1) ...
  +0.5*lambda1*((C(N+2,3:M)+C(N+1,3:M)).*U(N,2:M-1)...
  -(C(N+2,3:M)+2*C(N+1,3:M)+C(N,3:M)).*U(N,2:M-1)...
  +(C(N+1,3:M)+C(N,3:M)).*U(N-1,2:M-1)) ...
  +0.5*lambda2*((C(N+1,4:M+1)+C(N+1,3:M)).*U(N,3:M)...
  -(C(N,4:M+1)+2*C(N+1,3:M)+C(N,2:M-1)).*U(N,2:M-1)...
  +(C(N+1,2:M-1)+C(N+1,3:M)).*U(N,1:M-2));

% j=1 update;
  U_old(2:N-1, 1)=U(2:N-1, 1)...
  +0.5*lambda1*((C(4:N+1,2)+C(3:N,2)).*U(3:N,1)...
  -(C(4:N+1,2)+2*C(3:N,2)+C(2:N-1,2)).*U(2:N-1,1)...
  +(C(2:N-1,2)+C(3:N,2)).*U(1:N-2,1)) ...
  +0.5*lambda2*((C(3:N,3)+C(3:N,2)).*U(2:N-1,2)...
  -(C(3:N,3)+2*C(3:N,2)+C(3:N,1)).*U(2:N-1,1)...
  +(C(3:N,1)+C(3:N,2)).*U(2:N-1,1));

% j=M update
  U_old(2:N-1,M)=U(2:N-1,M)...
  +0.5*lambda1*((C(4:N+1,M+1)+C(3:N,M+1)).*U(3:N,M)...
  -(C(4:N+1,M+1)+2*C(3:N,M+1)+C(2:N-1,M+1)).*U(2:N-1,M)...
  +(C(2:N-1,M+1)+C(3:N,M+1)).*U(1:N-2,M))...
  +0.5*lambda2*((C(3:N,M+2)+C(3:N,M+1)).*U(2:N-1,M)...
  -(C(3:N,M+2)+2*C(3:N,M+1)+C(3:N,M)).*U(2:N-1,M)...
  +(C(3:N,M)+C(3:N,M+1)).*U(2:N-1,M-1));

 % i=1, j=1 update
  U_old(1, 1)=U(1, 1)...
  +0.5*lambda1*((C(3,2)+C(2,2)).*U(2,1)...
  -(C(3,2)+2*C(2,2)+C(1,2)).*U(1,1)...
  +(C(1,2)+C(2,2)).*U(1,1)) ...
  +0.5*lambda2*((C(2,3)+C(2,2)).*U(1,2)...
  -(C(2,3)+2*C(2,2)+C(2,1)).*U(1,1)...
  +(C(2,1)+C(2,2)).*U(1,1));

% i=N, j=1 update
  U_old(N, 1)=U(N, 1)...
  +0.5*lambda1*((C(N+2,2)+C(N+1,2)).*U(N,1)...
  -(C(N+2,2)+2*C(N+1,2)+C(N,2)).*U(N,1)...
  +(C(N,2)+C(N+1,2)).*U(N-1,1)) ...
  +0.5*lambda2*((C(N+1,3)+C(N+1,2)).*U(N,2)...
  -(C(N+1,3)+2*C(N+1,2)+C(N+1,1)).*U(N,1)...
  +(C(N+1,1)+C(N+1,2)).*U(N,1));

  

% j=M, i=1 update
  U_old(1, M)=U(1, M)...
  +0.5*lambda1*((C(3,M+1)+C(2,M+1)).*U(2,M)...
  -(C(3,M+1)+2*C(2,M+1)+C(1,M+1)).*U(1,M)...
  +(C(2,M+1)+C(1,M+1)).*U(1,M)) ...
  +0.5*lambda2*((C(2,M+2)+C(2,M+1)).*U(1,M)...
  -(C(2,M+2)+2*C(2,M+1)+C(2,M)).*U(1,M)...
  +(C(2,M)+C(2,M+1)).*U(1,M-1));

% j=M, i=N update
  U_old(N, M)=U(N, M) ...
  +0.5*lambda1*((C(N+2,M+1)+C(N+1,2)).*U(N,M)...
  -(C(N+2,M+1)+2*C(N+1,M+1)+C(N,M+1)).*U(N,M)...
  +(C(N,M+1)+C(N+1,M+1)).*U(N-1,M)) ...
  +0.5*lambda2*((C(N+1,M+2)+C(N+1,M+1)).*U(N,M)...
  -(C(N+1,M+2)+2*C(N+1,M+1)+C(N+1,M)).*U(N,M)...
  +(C(N+1,M)+C(N+1,M+1)).*U(N,M-1));

U=U_old;
 end
 
 %error 
 % exact solution
 
%  U_exact=zeros(N,M);
%  U_exact=exp(-2*pi^2*nu*dt*maxiter)*cos(pi.*X).*cos(pi.*Y);
%  U_exact=U_exact';
%  
%  norm(U_exact-U)*sqrt(dx*dy)
  %surf(U)
   snr(I2,U-I2)

subplot(2,2,1),imshow(I1),title('origional image');                              
subplot(2,2,2),imshow(J1),title('noisy image');
subplot(2,2,3),imshow(uint8(U)),title('denoise image');
subplot(2,2,4),imshow(uint8(U-I2)),title('noise');

% 