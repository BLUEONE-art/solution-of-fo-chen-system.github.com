clc;
clear;
%��ȡ����ԭͼһ�����ص��ļ���ͼ
I=imread('../author_images/Yun_Wu.jpg','jpg');         %��ȡͼ����Ϣ��lena
I=rgb2gray(I);
C=imresize(I,[375 300]);
imwrite(C,'../author_images/a4_Yun_Wu.png','png');
figure('color',[1 1 1]);imshow(C);title('�޸ĺ�ͼƬ');