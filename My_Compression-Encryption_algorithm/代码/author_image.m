clc;
clear;
%读取更改原图一个像素点后的加密图
I=imread('../author_images/Yun_Wu.jpg','jpg');         %读取图像信息，lena
I=rgb2gray(I);
C=imresize(I,[375 300]);
imwrite(C,'../author_images/a4_Yun_Wu.png','png');
figure('color',[1 1 1]);imshow(C);title('修改后图片');