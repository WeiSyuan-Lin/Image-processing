%add salt and papper noise,�¥��I�۶�
pic01=imread('f01.jpg');
H=rgb2gray(pic01)%���ܶ¥�mxn���x�}
H1=im2double(H) %�������double ���b0~1����
%imagesc(Graypic01)
R=rand(size(H1));
Z=H1+R;
B=find(Z<0.85);
W=find(Z>0.98);
for k=1:size(B)

Z(B(k))=0;
    
end

for k=1:size(W)

Z(W(k))=1;
    
end


%imshow(Z);
subplot(1,2,1);
imshow(pic01);

title('original photo');
hold on;

subplot(1,2,2);
imshow(Z);

title('modified photo');

hold off;
