pic01=imread('f01.jpg');
S=size(pic01);

pic01(round(S(1)/8):round(7*S(1)/8),round(1*S(2)/8):round(2*S(2)/8),1)=0
pic01(round(S(1)/8):round(7*S(1)/8),round(1*S(2)/8):round(2*S(2)/8),2)=0
pic01(round(S(1)/8):round(7*S(1)/8),round(1*S(2)/8):round(2*S(2)/8),3)=0


pic01(round(3*S(1)/8):round(4*S(1)/8),round(2*S(2)/8):round(6*S(2)/8),1)=0
pic01(round(3*S(1)/8):round(4*S(1)/8),round(2*S(2)/8):round(6*S(2)/8),2)=0
pic01(round(3*S(1)/8):round(4*S(1)/8),round(2*S(2)/8):round(6*S(2)/8),3)=0


pic01(round(S(1)/8):round(7*S(1)/8),round(6*S(2)/8):round(7*S(2)/8),1)=0
pic01(round(S(1)/8):round(7*S(1)/8),round(6*S(2)/8):round(7*S(2)/8),2)=0
pic01(round(S(1)/8):round(7*S(1)/8),round(6*S(2)/8):round(7*S(2)/8),3)=0


pic00=imread('f01.jpg');

subplot(1,2,1);
imshow(pic00);

title('original photo');
hold on;

subplot(1,2,2);
imshow(pic01);

title('modified photo');

hold off;