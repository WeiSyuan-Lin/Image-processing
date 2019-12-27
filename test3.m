
Im=imread('example.jpg');

I1=rgb2gray(Im);
I2=im2double(I1);
I=I2;


[N,M]=size(I);




nu=1;
dx=1/N;
dy=1/M;
dt=0.1*dx^2/nu;
delt=0.1;
lambda1=dt*nu/(dx)^2;
lambda2=dt*nu/(dy)^2;
% x=dx*(1/2:1:N-1/2);
% y=dy*(1/2:1:M-1/2);
C=zeros(N+2,M+2);
%coef=@(x) nu;
%coef=@(x) 1./sqrt(1+x);
K=3;
%coef=@(x) 1./(1+x./K);
coef=@(x)  exp(-x./K);

Ix=zeros(size(I));
Iy=zeros(size(I));
Ixx=zeros(size(I));
Iyy=zeros(size(I));
dLx=zeros(size(I));
dLy=zeros(size(I));
Ix1=zeros(size(I)); 
Iy1=zeros(size(I));
Ixx1=zeros(size(I));
Iyy1=zeros(size(I));


maxiter =500;
maxiter1=50;
maxiter2=2;

Inew=I;
mask=zeros(size(I));
mask(find(I==1))=1; %mask defined all black and inpainting region white

for itr = 1:maxiter
for k=1:maxiter1
    
 
Ix1(2:N-1,2:M-1)=(I(3:N,2:M-1)-I(1:N-2,2:M-1))/(2*dx);
Ix(find(mask==1))=Ix1(find(mask==1));

Iy1(2:N-1,2:M-1)=(I(2:N-1,3:M)-I(2:N-1,1:M-2))/(2*dy);
Iy(find(mask==1))=Iy1(find(mask==1));

 denN=sqrt((Ix.^2)+(Iy.^2)+0.000001);
 Nx=Ix./denN;
 Ny=Iy./denN;
 
    

 %Ixx(2:N-1,2:M-1)=(mask(3:N,2:M-1)+mask(1:N-2,2:M-1)-2*mask(2:N-1,2:M-1))/(dx.^2);
 %Iyy(2:N-1,2:M-1)=(mask(2:N-1,3:M)+mask(2:N-1,1:M-2)-2*mask(2:N-1,2:M-1))/(dy.^2);
%  
Ixx1(2:N-1,2:M-1)=(I(3:N,2:M-1)+I(1:N-2,2:M-1)-2*I(2:N-1,2:M-1))/(dx.^2);
Ixx(find(mask==1))=Ixx1(find(mask==1)); 
 
 Iyy1(2:N-1,2:M-1)=(I(2:N-1,3:M)+I(2:N-1,1:M-2)-2*I(2:N-1,2:M-1))/(dx.^2);
Iyy(find(mask==1))=Iyy1(find(mask==1));  

 L=(Ixx.^2+Iyy.^2);
 dLx(2:N-1,2:M-1)=(dx.^2)*(dy.^2)*(L(3:N,2:M-1)-L(1:N-2,2:M-1));
 dLy(2:N-1,2:M-1)=(dy.^2)*(dx.^2)*(L(2:N-1,3:M)-L(2:N-1,1:M-2));
 
 B=0.2*((dLx.*(-Ny))+(dLy.*Nx));
 Ixb=zeros(size(I));
 Iyb=zeros(size(I));
 Ixf=zeros(size(I));
 Iyf=zeros(size(I));
 
 
 Ixb(2:N,2:M)=(I(2:N,2:M)-I(1:N-1,2:M))/dx;         
 Iyb(2:N,2:M)=(I(2:N,2:M)-I(2:N,1:M-1))/dy;  
 Ixbm=min(Ixb,0);
 Iybm=min(Iyb,0);
 IxbM=max(Ixb,0);
 IybM=max(Iyb,0);

 
 Ixf(1:N-1,1:M-1)=(I(2:N,1:M-1)-I(1:N-1,1:M-1))/dx;         
 Iyf(1:N-1,1:M-1)=(I(1:N-1,2:M)-I(1:N-1,1:M-1))/dy;  
 Ixfm=min(Ixf,0);
 Iyfm=min(Iyf,0);
 IxfM=max(Ixf,0);
 IyfM=max(Iyf,0);

 G1=sqrt(Ixbm.^2+IxfM.^2+Iybm.^2+IyfM.^2);
 G2=sqrt(IxbM.^2+Ixfm.^2+IybM.^2+Iyfm.^2);
 G=G2;
 G(find(B>0))=G1(find(B>0));
 It=B.*G;
 Inew=I+delt*It;
 I=Inew;
%  
  I(find(I>1))=1; 
  I(find(I<0))=0;

end
U=I;  

for k=1:maxiter2
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

        U_old(2:N-1, 2:M-1)=U_old(2:N-1, 2:M-1);

        % i=1 update
        U_old(1, 2:M-1)=U(1, 2:M-1) ...
          +0.5*lambda1*((C(3,3:M)+C(2,3:M)).*U(2,2:M-1)...
          -(2*C(2,3:M)+C(3,3:M)+C(1,3:M)).*U(1,2:M-1)...
                      +(C(1,3:M)+C(2,3:M)).*U(1,2:M-1)) ...
          +0.5*lambda2*((C(2,4:M+1)+C(2,3:M)).*U(1,3:M)...
          -(C(2,4:M+1)+2*C(2,3:M)+C(2,2:M-1)).*U(1,2:M-1)...
            +(C(2,2:M-1)+C(2,3:M)).*U(1,1:M-2));

        U_old(1, 2:M-1)=U_old(1, 2:M-1) ;

        % i=N update
         U_old(N, 2:M-1)=U(N, 2:M-1) ...
          +0.5*lambda1*((C(N+2,3:M)+C(N+1,3:M)).*U(N,2:M-1)...
          -(2*C(N+1,3:M)+C(N+2,3:M)+C(N,3:M)).*U(N,2:M-1)...
                      +(C(N,3:M)+C(N+1,3:M)).*U(N-1,2:M-1)) ...
          +0.5*lambda2*((C(N+1,4:M+1)+C(N+1,3:M)).*U(N,3:M)...
          -(C(N+1,4:M+1)+2*C(N+1,3:M)+C(N+1,2:M-1)).*U(N,2:M-1)...
            +(C(N+1,2:M-1)+C(N+1,3:M)).*U(N,1:M-2));

        U_old(N, 2:M-1)=U_old(N, 2:M-1) ;

        % j=1 update
         U_old(2:N-1, 1)=U(2:N-1, 1) ...
          +0.5*lambda1*((C(4:N+1,2)+C(3:N,2)).*U(3:N,1)...
          -(2*C(3:N,2)+C(4:N+1,2)+C(2:N-1,2)).*U(2:N-1,1)...
                      +(C(2:N-1,2)+C(3:N,2)).*U(1:N-2,1)) ...
          +0.5*lambda2*((C(3:N,3)+C(3:N,2)).*U(2:N-1,2)...
          -(C(3:N,3)+2*C(3:N,2)+C(3:N,1)).*U(2:N-1,1)...
            +(C(3:N,1)+C(3:N,2)).*U(2:N-1,1));

         U_old(2:N-1, 1)=U_old(2:N-1, 1) ;


        % j=M update
         U_old(2:N-1, M)=U(2:N-1, M) ...
          +0.5*lambda1*((C(4:N+1,M+1)+C(3:N,M+1)).*U(3:N,M)...
          -(2*C(3:N,M+1)+C(4:N+1,M+1)+C(2:N-1,M+1)).*U(2:N-1,M)...
                      +(C(2:N-1,M+1)+C(3:N,M+1)).*U(1:N-2,M)) ...
          +0.5*lambda2*((C(3:N,M+2)+C(3:N,M+1)).*U(2:N-1,M)...
          -(C(3:N,M+2)+2*C(3:N,M+1)+C(3:N,M)).*U(2:N-1,M)...
            +(C(3:N,M)+C(3:N,M+1)).*U(2:N-1,M-1));

        U_old(2:N-1, M)=U_old(2:N-1, M) ;


         % i=1, j=1 update
        U_old(1, 1)=U(1, 1) ...
          +0.5*lambda1*((C(3,2)+C(2,2)).*U(2,1)...
          -(2*C(2,2)+C(3,2)+C(1,2)).*U(1,1)...
                      +(C(1,2)+C(2,2)).*U(1,1)) ...
          +0.5*lambda2*((C(2,3)+C(2,2)).*U(1,2)...
          -(C(2,3)+2*C(2,2)+C(2,1)).*U(1,1)...
            +(C(2,2)+C(2,1)).*U(1,1));

        U_old(1, 1)=U_old(1, 1) ;

        % i=N, j=1 update
         U_old(N, 1)=U(N, 1) ...
          +0.5*lambda1*((C(N+2,2)+C(N+1,2)).*U(N,1)...
          -(2*C(N+1,2)+C(N+2,2)+C(N,2)).*U(N,1)...
                      +(C(N,2)+C(N+1,2)).*U(N-1,1)) ...
          +0.5*lambda2*((C(N+1,3)+C(N+1,2)).*U(N,2)...
          -(C(N+1,3)+2*C(N+1,2)+C(N+1,1)).*U(N,1)...
            +(C(N+1,2)+C(N+1,1)).*U(N,1));

        U_old(N, 1)=U_old(N, 1) ;  

        % j=M, i=1 update
         U_old(1, M)=U(1, M) ...
          +0.5*lambda1*((C(3,M+1)+C(2,M+1)).*U(2,M)...
          -(2*C(2,M+1)+C(3,M+1)+C(1,M+1)).*U(1,M)...
                      +(C(1,M+1)+C(2,M+1)).*U(1,M)) ...
          +0.5*lambda2*((C(2,M+2)+C(2,M+1)).*U(1,M)...
          -(C(2,M+2)+2*C(2,M+1)+C(2,M)).*U(1,M)...
            +(C(2,M+1)+C(2,M)).*U(1,M-1));
        U_old(1, M)=U_old(1, M) ;

        % j=M, i=N update
         U_old(N, M)=U(N, M) ...
          +0.5*lambda1*((C(N+2,M+1)+C(N+1,M+1)).*U(N,M)...
          -(2*C(N+1,M+1)+C(N+2,M+1)+C(N,M+1)).*U(N,M)...
                      +(C(N,M+1)+C(N+1,M+1)).*U(N-1,M)) ...
          +0.5*lambda2*((C(N+1,M+2)+C(N+1,M+1)).*U(N,M)...
          -(C(N+1,M+2)+2*C(N+1,M+1)+C(N+1,M)).*U(N,M)...
            +(C(N+1,M+1)+C(N+1,M)).*U(N,M-1));
        U_old(N, M)=U_old(N, M);

        U=U_old;
        
end
end
% 

 for i=1:N
       for j=1:M
           if (mask(i,j)==1)
               if U(i,j)>1
                   U(i,j)=1;
               end
               if U(i,j)<0
                    U(i,j)=0;
               end
               I(i,j) = U(i,j);
           end
       end
 end
 
 imshow(I);
