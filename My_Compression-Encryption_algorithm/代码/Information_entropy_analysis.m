%% 加密系统信息熵分析
%   @author:董昊
%   @date:2020.04.17
%% 原始图片R,G,B通道信息熵
clear;clc;
% I=imread('../原始、加密、解密图片/lena/lena.png','png');         %读取图像信息
% I=imread('../原始、加密、解密图片/辣椒/peppers.png','png');       %读取图像信息，辣椒
% I=imread('../原始、加密、解密图片/狒狒/baboon.png','png');         %读取图像信息，狒狒
% I=imread('../原始、加密、解密图片/飞机/airplane.png','png');         %读取图像信息，飞机
% I=imread('../原始、加密、解密图片/4.1.05房子/4.1.05.tiff','tiff');         %读取图像信息
I=imread('../原始、加密、解密图片/4.2.06小船/4.2.06.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/汽车/house.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/城堡/kodim08.png','png');         %读取图像信息

% J1=imread('../原始、加密、解密图片/lena/加密后的lena.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/辣椒/加密后的peppers.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/狒狒/加密后的baboon.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/飞机/加密后的airplane.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/4.1.05房子/加密后的4.1.05.png','png');             %读取图像信息
J1=imread('../原始、加密、解密图片/4.2.06小船/加密后的4.2.06.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/汽车/加密后的house.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/城堡/加密后的kodim08.png','png');             %读取图像信息

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
Q_R=J1(:,:,1);
Q_G=J1(:,:,2);
Q_B=J1(:,:,3);
%R通道
I1=uint8(I1);
I2=uint8(I2);
I3=uint8(I3);
T1_R=imhist(I1);   %统计图像R（Y）通道灰度值从0~255的分布情况，存至T1
S1_R=sum(T1_R);     %计算整幅图像R通道的灰度值
xxs1_R=0;           %初始化原始图片R通道相关性为0（信息熵）
%G通道
T1_G=imhist(I2);
S1_G=sum(T1_G);
xxs1_G=0;
%B通道
T1_B=imhist(I3);
S1_B=sum(T1_B);
xxs1_B=0;

for i=1:256
    pp1_R=T1_R(i)/S1_R;   %每个灰度值占比，即每个灰度值的概率
    pp1_G=T1_G(i)/S1_G;
    pp1_B=T1_B(i)/S1_B;
    if pp1_R~=0
        xxs1_R=xxs1_R-pp1_R*log2(pp1_R);%计算原始图像R通道的相关性（R通道信息熵）
    end
    if pp1_G~=0
        xxs1_G=xxs1_G-pp1_G*log2(pp1_G);%计算原始图像G通道的相关性（G通道信息熵）
    end
    if pp1_B~=0
        xxs1_B=xxs1_B-pp1_B*log2(pp1_B);%计算原始图像B通道的相关性（B通道信息熵）
    end
end
%% 加密后信息熵
Q_R1=uint8(Q_R);
Q_G1=uint8(Q_G);
Q_B1=uint8(Q_B);
%R通道
T2_R=imhist(Q_R1);
S2_R=sum(T2_R);
xxs2_R=0;
%G通道
T2_G=imhist(Q_G1);
S2_G=sum(T2_G);
xxs2_G=0;
%B通道
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
disp('信息熵：');
disp(['原始图片R通道信息熵=',num2str(xxs1_R),'  原始图片G通道信息熵=',num2str(xxs1_G),'  原始图片B通道信息熵=',num2str(xxs1_B)]);
disp(['加密图片R通道信息熵=',num2str(xxs2_R),'  加密图片G通道信息熵=',num2str(xxs2_G),'  加密图片B通道信息熵=',num2str(xxs2_B)]);