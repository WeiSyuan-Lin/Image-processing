% 2D heat solver 

clc
clear all;
% numerical parameters
% physcial parameters

I=imread('f01.jpg');
I1=rgb2gray(I);
I2=im2double(I1);

J1=imnoise(I2,'salt & pepper',0.02);

J2=double(J1);

nu=1;
[N,M]=size(I2);

dx=1/N;
dy=1/M;
dt=0.1*dx^2/nu;
lambda1=dt*nu/(dx)^2;
lambda2=dt*nu/(dy)^2;
maxiter=100;%1/dt;
x=dx*(1/2:1:N-1/2);
y=dy*(1/2:1:M-1/2);

% initial data

U=zeros(N,M);
U=J2;

 for k=1:maxiter
 U_old=U;
 
 U_old(2:N-1, 2:M-1)=U(2:N-1, 2:M-1)...
  +lambda1*(U(3:N,2:M-1)-2*U(2:N-1,2:M-1)+U(1:N-2,2:M-1)) ...
  +lambda2*(U(2:N-1,3:M)-2*U(2:N-1,2:M-1)+U(2:N-1,1:M-2));

% i=1 update
U_old(1, 2:M-1)=U(1, 2:M-1) ...
  +lambda1*(U(2,2:M-1)-2*U(1,2:M-1)+U(1,2:M-1)) ...
  +lambda2*(U(1,3:M)-2*U(1,2:M-1)+U(1,1:M-2));

% i=N update
 U_old(N, 2:M-1)=U(N, 2:M-1) ...
  +lambda1*(U(N,2:M-1)-2*U(N,2:M-1)+U(N-1,2:M-1)) ...
  +lambda2*(U(N,3:M)-2*U(N,2:M-1)+U(N,1:M-2));
 
% j=1 update
 U_old(2:N-1, 1)=U(2:N-1, 1) ...
  +lambda1*(U(3:N,1)-2*U(2:N-1,1)+U(1:N-2,1)) ...
  +lambda2*(U(2:N-1,2)-2*U(2:N-1,1)+U(2:N-1,1));    

% j=M update
 U_old(2:N-1, M)=U(2:N-1, M) ...
  +lambda1*(U(3:N,M)-2*U(2:N-1,M)+U(1:N-2,M)) ...
  +lambda2*(U(2:N-1,M)-2*U(2:N-1,M)+U(2:N-1,M-1));

 % i=1, j=1 update
U_old(1, 1)=U(1, 1) ...
  +lambda1*(U(2,1)-2*U(1,1)+U(1,1)) ...
  +lambda2*(U(1,2)-2*U(1,1)+U(1,1));
% i=N, j=1 update
 U_old(N, 1)=U(N, 1) ...
  +lambda1*(U(N,1)-2*U(N,1)+U(N-1,1)) ...
  +lambda2*(U(N,2)-2*U(N,1)+U(N,1));
 
  

% j=M, i=1 update
 U_old(1, M)=U(1, M) ...
  +lambda1*(U(2,M)-2*U(1,M)+U(1,M)) ...
  +lambda2*(U(1,M)-2*U(1,M)+U(1,M-1));

% j=M, i=N update
 U_old(N, M)=U(N, M) ...
  +lambda1*(U(N,M)-2*U(N,M)+U(N-1,M)) ...
  +lambda2*(U(N,M)-2*U(N,M)+U(N,M-1));

U=U_old;
 end
 
 %error 
 snr(I2,U-I2)

subplot(2,2,1),imshow(I1),title('origional image');                              
subplot(2,2,2),imshow(J1),title('noisy image');
subplot(2,2,3),imshow(uint8(U)),title('denoise image');
subplot(2,2,4),imshow(uint8(U-I2)),title('noise');
% 