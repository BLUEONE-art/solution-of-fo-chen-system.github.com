%% 加密系统NPCR（Number of pixels change rate）&UACI（unifed average changing intensity）分析
%   @author:董昊
%   @date:2020.04.17
clc;
clear;
%读取更改原图一个像素点后的加密图
image1=imread('../原始、加密、解密图片/lena/改变一个像素加密后的lena.png','png');             %读取图像信息
% image1=imread('../原始、加密、解密图片/辣椒/改变一个像素加密后的peppers.png','png');             %读取图像信息
% image1=imread('../原始、加密、解密图片/狒狒/改变一个像素加密后的baboon.png','png');             %读取图像信息
% image1=imread('../原始、加密、解密图片/飞机/改变一个像素加密后的airplane.png','png');             %读取图像信息
% image1=imread('../原始、加密、解密图片/汽车/改变一个像素加密后的house.png','png');             %读取图像信息
% image1=imread('../原始、加密、解密图片/4.2.06小船/改变一个像素加密后的4.2.06.png','png');             %读取图像信息

I1=image1(:,:,1);
I2=image1(:,:,2);
I3=image1(:,:,3);
%读取没有更改原图的加密图
image2=imread('../原始、加密、解密图片/lena/加密后的lena.png','png');              %读取图像信息
% image2=imread('../原始、加密、解密图片/辣椒/加密后的peppers.png','png');              %读取图像信息
% image2=imread('../原始、加密、解密图片/狒狒/加密后的baboon.png','png');              %读取图像信息
% image2=imread('../原始、加密、解密图片/飞机/加密后的airplane.png','png');              %读取图像信息
% image2=imread('../原始、加密、解密图片/汽车/加密后的house.png','png');              %读取图像信息
% image2=imread('../原始、加密、解密图片/4.2.06小船/加密后的4.2.06.png','png');              %读取图像信息

J1=image2(:,:,1);
J2=image2(:,:,2);
J3=image2(:,:,3);
%显示图片
figure;imshow(image2),title('原图像');
figure;imshow(image1),title('改变后图像');
%计算相同位置灰度值相等的个数
[M,N]=size(I1);
m=0;
n=0;
p=0;
u1=0;
u2=0;
u3=0;
for i=1:M
    for j=1:N
        u1=u1+abs(double(I1(i,j))-double(J1(i,j)))/255;
        u2=u2+abs(double(I2(i,j))-double(J2(i,j)))/255;
        u3=u3+abs(double(I3(i,j))-double(J3(i,j)))/255;
        if I1(i,j)~=J1(i,j)
            m=m+1;
        end
        if I2(i,j)~=J2(i,j)
            n=n+1;
        end
        if I3(i,j)~=J3(i,j)
            p=p+1;
        end
    end
end
ua=sum(abs(double(I1(:))-double(J1(:)))/255);%简单的式子，一个就搞定
%计算NPCR
NPCR1=(m*100)/(M*N);
NPCR2=(n*100)/(M*N);
NPCR3=(p*100)/(M*N);
NPCR=sum(double(I1(:))~=double(J1(:)))*100/(M*N);%简单的方法
UACI1=(u1*100)/(M*N);
UACI2=(u2*100)/(M*N);
UACI3=(u3*100)/(M*N);
UACI=(ua*100)/(M*N);
disp(['R通道的NPCR值为：',num2str(NPCR1)]);
disp(['G通道的NPCR值为：',num2str(NPCR2)]);
disp(['B通道的NPCR值为：',num2str(NPCR3)]);
% disp(['NPCR值为：',num2str(NPCR1)]);
disp(['R通道的UACI值为：',num2str(UACI1)]);%方法一
disp(['G通道的UACI值为：',num2str(UACI2)]);%方法一
disp(['B通道的UACI值为：',num2str(UACI3)]);%方法一
% disp(['UACI值为：',num2str(UACI1)]);%方法二