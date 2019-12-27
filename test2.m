Im=imread('example.jpg');

I1=rgb2gray(Im);
I2=im2double(I1);
I=I2;


[N,M]=size(I);

dx=1/N;
dy=1/M;
%dx=1;
%dy=1;
dt=0.1;
maxiter=20;


%imshow((mask));
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
Inew=I;
%I=imgaussfilt(I);
%for n=1:5
mask=zeros(size(I));
mask(find(I==1))=1; %mask defined all black and inpainting region white

for k=1:maxiter
    
   
%外部顏色進不去，要用原圖在mask附近的值
% Ix(2:N-1,2:M-1)=(mask(3:N,2:M-1)-mask(1:N-2,2:M-1))/dx;
% Iy(2:N-1,2:M-1)=(mask(2:N-1,3:M)-mask(2:N-1,1:M-2))/dy;         
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
 Inew=I+dt*It;
 I=Inew;
%  
  I(find(I>1))=1; 
  I(find(I<0))=0;
end


%end
I=imgaussfilt(I);
imshow(I);