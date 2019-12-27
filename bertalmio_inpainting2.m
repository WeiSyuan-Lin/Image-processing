clear all;clc;
img=[128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ;128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 ; 128 128 128 128 128 128 128 128 128 128 255 255 255 255 255 128 128 128 128 128 128 128 128 128 128; 128 128 128 128 128 128 128 128 128 128 255 255 255 255 255 128 128 128 128 128 128 128 128 128 128; 128 128 128 128 128 128 128 128 128 128 255 255 255 255 255 128 128 128 128 128 128 128 128 128 128; 0 0 0 0 0 0 0 0 0 0 255 255 255 255 255 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 255 255 255 255 255 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Im=imread('example.jpg');

I1=rgb2gray(Im);
I2=double(I1);
img=I2;
mask=zeros(size(img));
mask(img==255)=255; %mask defined all black and inpainting region white
% possible start of iteration
for k=1:200
% bound=edge(mask,'log'); % detect boundary
for x=3:(size(img,1)-2)
    for y=3:(size(img,2)-2)
        if (mask(x,y)==255) % if it belongs to the mask
%             mask(x-1,y)=0;mask(x+1,y)=0;mask(x,y-1)=0;mask(x,y+1)=0;mask(x,y)=0;
            %3.4
            Ix=(img(x+1,y)-img(x-1,y))/2;
            Iy=(img(x,y+1)-img(x,y-1))/2;
            denN=sqrt((Ix*Ix)+(Iy*Iy)+0.000001); %epsilon 0.000001 to aviod divide by zero
            Nx=Ix/denN;
            Ny=Iy/denN;
            %3.3
                 p1=x+1;q1=y;
                 Lnew1=img(p1+1,q1)+img(p1-1,q1)+img(p1,q1+1)+img(p1,q1-1)-4*img(p1,q1);
                 p2=x-1;q2=y;
                 Lnew2=img(p2+1,q2)+img(p2-1,q2)+img(p2,q2+1)+img(p2,q2-1)-4*img(p2,q2);
                 p3=x;q3=y+1;
                 Lnew3=img(p3+1,q3)+img(p3-1,q3)+img(p3,q3+1)+img(p3,q3-1)-4*img(p3,q3);
                 p4=x;q4=y-1;
                 Lnew4=img(p4+1,q4)+img(p4-1,q4)+img(p4,q4+1)+img(p4,q4-1)-4*img(p4,q4);
                 deLx=Lnew1-Lnew2;
                 deLy=Lnew3-Lnew4;
            %3.6
            beta=0.05*((deLx*(-Ny))+(deLy*Nx));
            %3.5
            Ixf=img(x+1,y)-img(x,y);
            Ixb=img(x,y)-img(x-1,y);
            Iyf=img(x,y+1)-img(x,y);
            Iyb=img(x,y)-img(x,y-1);
            Ixfm=min(Ixf,0);
            IxfM=max(Ixf,0);
            Ixbm=min(Ixb,0);
            IxbM=max(Ixb,0);
            Iyfm=min(Iyf,0);
            IyfM=max(Iyf,0);
            Iybm=min(Iyb,0);
            IybM=max(Iyb,0);
            if(beta>0)
                ModeLI=sqrt((Ixbm*Ixbm)+(IxfM*IxfM)+(Iybm*Iybm)+(IyfM*IyfM));
            elseif (beta<=0)
                ModeLI=sqrt((IxbM*IxbM)+(Ixfm*Ixfm)+(IybM*IybM)+(Iyfm*Iyfm));
            end
            %3.2
            It=beta*ModeLI;
            %3.1
            img(x,y)=img(x,y)+0.1*(It);
            img=double(uint8(int64(img)));% the values of It goes in e+05 and increases at every iteration so this statement to keep it to 255

        end
    end
end
end

imshow(uint8(img))