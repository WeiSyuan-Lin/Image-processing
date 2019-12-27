Im=imread('example.jpg');

I1=rgb2gray(Im);
I2=double(I1);
I=I2;


[N,M]=size(I);

dx=1/N;
dy=1/M;
dt=0.1;
maxiter=20;


mask=zeros(size(I));
%mask(I==255)=255; %mask defined all black and inpainting region white
mask(find(245<I&I<256))=255;
for k=1:maxiter
    
    for i=2:N-1
        for j=2:M-1
    
    if (mask(i,j)==255)

        Ix=(I(i+1,j)-I(i-1,j))/2; 
        Iy=(I(i,j+1)-I(i,j-1))/2;
        denN=sqrt((Ix*Ix)+(Iy*Iy)+0.000001); %epsilon 0.000001 to aviod divide by zero
            Nx=Ix/denN;
            Ny=Iy/denN;



                 p1=i+1;q1=j;
                 L1=(I(p1+1,q1)+I(p1-1,q1)-2*I(p1,q1))+(I(p1,q1+1)+I(p1,q1-1)-2*I(p1,q1));
                 p2=i-1;q2=j;
                 L2=(I(p2+1,q2)+I(p2-1,q2)-2*I(p2,q2))+(I(p2,q2+1)+I(p2,q2-1)-2*I(p2,q2));
                 p3=i;q3=j+1;
                 L3=(I(p3+1,q3)+I(p3-1,q3)-2*I(p3,q3))+(I(p3,q3+1)+I(p3,q3-1)-2*I(p3,q3));
                 p4=i;q4=j-1;
                 L4=(I(p4+1,q4)+I(p4-1,q4)-2*I(p4,q4))+(I(p4,q4+1)+I(p4,q4-1)-2*I(p4,q4));
                 

%L(i,j)=(Ixx(i,j).^2+Iyy(i,j).^2)*((dx^2)*(dy^2));
dLx=(L1-L2);
dLy=(L3-L4);


% LN(i,j)=sqrt(Ix(i,j).^2+Iy(i,j).^2+0.00001);
%  B(i,j)=(Ix(i,j).*dLy(i,j)-Iy(i,j).*dLx(i,j))./LN(i,j);
 B=0.05*((dLx*(-Ny))+(dLy*Nx));
            Ixf=(I(i+1,j)-I(i,j));
            Ixb=(I(i,j)-I(i-1,j));
            Iyf=(I(i,j+1)-I(i,j));
            Iyb=(I(i,j)-I(i,j-1));
            
            Ixfm=min(Ixf,0);
            IxfM=max(Ixf,0);
            Ixbm=min(Ixb,0);
            IxbM=max(Ixb,0);
            Iyfm=min(Iyf,0);
            IyfM=max(Iyf,0);
            Iybm=min(Iyb,0);
            IybM=max(Iyb,0);
           
            
              if(B>0)
                ModeLI=sqrt((Ixbm*Ixbm)+(IxfM*IxfM)+(Iybm*Iybm)+(IyfM*IyfM));
            elseif (B<=0)
                ModeLI=sqrt((IxbM*IxbM)+(Ixfm*Ixfm)+(IybM*IybM)+(Iyfm*Iyfm));
            end
          It=B*ModeLI;
      I(i,j)=I(i,j)+dt*It;  
     I=double(uint8(int64(I)));
    end
        end
        
    end
    
end
imshow(uint8(I));
