%% ����ϵͳ��Ϣ�ط���
%   @author:���
%   @date:2020.04.17
%% ԭʼͼƬR,G,Bͨ����Ϣ��
clear;clc;
% I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/peppers.png','png');       %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/baboon.png','png');         %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/airplane.png','png');         %��ȡͼ����Ϣ���ɻ�
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/4.1.05.tiff','tiff');         %��ȡͼ����Ϣ
I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/4.2.06.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/house.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/kodim08.png','png');         %��ȡͼ����Ϣ

% J1=imread('../ԭʼ�����ܡ�����ͼƬ/lena/���ܺ��lena.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��peppers.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��baboon.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/���ܺ��airplane.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/���ܺ��4.1.05.png','png');             %��ȡͼ����Ϣ
J1=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/���ܺ��4.2.06.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��house.png','png');             %��ȡͼ����Ϣ
% J1=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/���ܺ��kodim08.png','png');             %��ȡͼ����Ϣ

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
Q_R=J1(:,:,1);
Q_G=J1(:,:,2);
Q_B=J1(:,:,3);
%Rͨ��
I1=uint8(I1);
I2=uint8(I2);
I3=uint8(I3);
T1_R=imhist(I1);   %ͳ��ͼ��R��Y��ͨ���Ҷ�ֵ��0~255�ķֲ����������T1
S1_R=sum(T1_R);     %��������ͼ��Rͨ���ĻҶ�ֵ
xxs1_R=0;           %��ʼ��ԭʼͼƬRͨ�������Ϊ0����Ϣ�أ�
%Gͨ��
T1_G=imhist(I2);
S1_G=sum(T1_G);
xxs1_G=0;
%Bͨ��
T1_B=imhist(I3);
S1_B=sum(T1_B);
xxs1_B=0;

for i=1:256
    pp1_R=T1_R(i)/S1_R;   %ÿ���Ҷ�ֵռ�ȣ���ÿ���Ҷ�ֵ�ĸ���
    pp1_G=T1_G(i)/S1_G;
    pp1_B=T1_B(i)/S1_B;
    if pp1_R~=0
        xxs1_R=xxs1_R-pp1_R*log2(pp1_R);%����ԭʼͼ��Rͨ��������ԣ�Rͨ����Ϣ�أ�
    end
    if pp1_G~=0
        xxs1_G=xxs1_G-pp1_G*log2(pp1_G);%����ԭʼͼ��Gͨ��������ԣ�Gͨ����Ϣ�أ�
    end
    if pp1_B~=0
        xxs1_B=xxs1_B-pp1_B*log2(pp1_B);%����ԭʼͼ��Bͨ��������ԣ�Bͨ����Ϣ�أ�
    end
end
%% ���ܺ���Ϣ��
Q_R1=uint8(Q_R);
Q_G1=uint8(Q_G);
Q_B1=uint8(Q_B);
%Rͨ��
T2_R=imhist(Q_R1);
S2_R=sum(T2_R);
xxs2_R=0;
%Gͨ��
T2_G=imhist(Q_G1);
S2_G=sum(T2_G);
xxs2_G=0;
%Bͨ��
T2_B=imhist(Q_B1);
S2_B=sum(T2_B);
xxs2_B=0;
for i=1:256
    pp2_R=T2_R(i)/S2_R;
    pp2_G=T2_G(i)/S2_R;
    pp2_B=T2_B(i)/S2_B;
    if pp2_R~=0
        xxs2_R=xxs2_R-pp2_R*log2(pp2_R);
    end
    if pp2_G~=0
        xxs2_G=xxs2_G-pp2_G*log2(pp2_G);
    end
    if pp2_B~=0
        xxs2_B=xxs2_B-pp2_B*log2(pp2_B);
    end
end
disp('��Ϣ�أ�');
disp(['ԭʼͼƬRͨ����Ϣ��=',num2str(xxs1_R),'  ԭʼͼƬGͨ����Ϣ��=',num2str(xxs1_G),'  ԭʼͼƬBͨ����Ϣ��=',num2str(xxs1_B)]);
disp(['����ͼƬRͨ����Ϣ��=',num2str(xxs2_R),'  ����ͼƬGͨ����Ϣ��=',num2str(xxs2_G),'  ����ͼƬBͨ����Ϣ��=',num2str(xxs2_B)]);