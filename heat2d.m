% 2D heat solver 

clc
clear all;
% numerical parameters
% physcial parameters
nu=1;
N=40;
M=40;
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





% construction of the matrix
% d1=[1-lambda1-lambda2 repmat(1-2*lambda1-lambda2,[1 N-2]) 1-lambda1-lambda2];
% d2=[1-lambda1-2*lambda2 repmat(1-2*lambda1-2*lambda2,[1 N-2]) 1-lambda1-2*lambda2]; 
% 
% A=diag([d1 repmat(d2, [1 M-2]) d1])+lambda1*diag(ones(N*M-1,1),1)...
%                                   +lambda1*diag(ones(N*M-1,1),-1)...
%                                   +lambda2*diag(ones(N*M-N,1),N)...
%                                   +lambda2*diag(ones(N*M-N,1),-N);
                              
                              
%  A2=lambda1*diag(ones(N-1,1),1)+diag([-lambda1 repmat(-2*lambda1,[1 N-2]) -lambda1]) ...
%      +lambda1*diag(ones(N-1,1),-1);
%  A3=lambda2*diag(ones(M-1,1),1)+diag([-lambda2 repmat(-2*lambda2,[1 M-2]) -lambda2]) ...
%      +lambda2*diag(ones(M-1,1),-1);
%  A4=eye(N*M)+kron(eye(M), A2) + kron(A3, eye(N));
 
 for k=1:maxiter
 U_old=U;
 
 U_old(2:N-1, 2:M-1)=U(2:N-1, 2:M-1) ...
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
 
 % exact solution
 
 U_exact=zeros(N,M);
 U_exact=exp(-2*pi^2*nu*dt*maxiter)*cos(pi.*X).*cos(pi.*Y);
 U_exact=U_exact';
 
 norm(U_exact-U)*sqrt(dx*dy)
 surf(U)
                              
                              


% 