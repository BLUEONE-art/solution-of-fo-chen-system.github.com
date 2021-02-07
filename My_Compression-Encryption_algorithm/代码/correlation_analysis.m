%% 加密系统像素相关性分析
%   @author:董昊
%   @date:2020.04.17
%% 原始图像相邻像素相关性分析
%{
先随机在0~M-1行和0~N-1列选中5000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
clear;clc;
I=imread('../原始、加密、解密图片/lena/lena.png','png');         %读取图像信息
% I=imread('../原始、加密、解密图片/辣椒/peppers.png','png');       %读取图像信息，辣椒
% I=imread('../原始、加密、解密图片/狒狒/baboon.png','png');         %读取图像信息，狒狒
% I=imread('../原始、加密、解密图片/飞机/airplane.png','png');         %读取图像信息，飞机
% I=imread('../原始、加密、解密图片/4.1.05房子/4.1.05.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/4.2.06小船/4.2.06.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/汽车/house.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/城堡/kodim08.png','png');         %读取图像信息

J1=imread('../原始、加密、解密图片/lena/加密后的lena.png','png');             %读取图像信息
% J1=imread('../原始、加密、解密图片/辣椒/加密后的peppers.png','png');       %读取图像信息，辣椒
% J1=imread('../原始、加密、解密图片/狒狒/加密后的baboon.png','png');         %读取图像信息，狒狒
% J1=imread('../原始、加密、解密图片/飞机/加密后的airplane.png','png');         %读取图像信息，飞机
% J1=imread('../原始、加密、解密图片/4.1.05房子/加密后的4.1.05.png','png');         %读取图像信息
% J1=imread('../原始、加密、解密图片/4.2.06小船/加密后的4.2.06.png','png');         %读取图像信息
% J1=imread('../原始、加密、解密图片/汽车/加密后的house.png','png');         %读取图像信息
% J1=imread('../原始、加密、解密图片/城堡/加密后的kodim08.png','png');         %读取图像信息

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
Q_R=J1(:,:,1);
Q_G=J1(:,:,2);
Q_B=J1(:,:,3);
[M,N]=size(I1);                      %将图像的行列赋值给M,N
NN=10000;    %随机取5000对像素点
x1=ceil(rand(1,NN)*(M-1));      %生成5000个1~M-1的随机整数作为行
y1=ceil(rand(1,NN)*(N-1));      %生成5000个1~N-1的随机整数作为列
%预分配内存
XX_R_SP=zeros(1,NN);YY_R_SP=zeros(1,NN);        %水平
XX_G_SP=zeros(1,NN);YY_G_SP=zeros(1,NN);
XX_B_SP=zeros(1,NN);YY_B_SP=zeros(1,NN);
XX_R_CZ=zeros(1,NN);YY_R_CZ=zeros(1,NN);        %垂直
XX_G_CZ=zeros(1,NN);YY_G_CZ=zeros(1,NN);
XX_B_CZ=zeros(1,NN);YY_B_CZ=zeros(1,NN);
XX_R_DJX=zeros(1,NN);YY_R_DJX=zeros(1,NN);      %对角线
XX_G_DJX=zeros(1,NN);YY_G_DJX=zeros(1,NN);
XX_B_DJX=zeros(1,NN);YY_B_DJX=zeros(1,NN);
for i=1:NN
    %水平
    XX_R_SP(i)=I1(x1(i),y1(i));
    YY_R_SP(i)=I1(x1(i)+1,y1(i));
    XX_G_SP(i)=I2(x1(i),y1(i));
    YY_G_SP(i)=I2(x1(i)+1,y1(i));
    XX_B_SP(i)=I3(x1(i),y1(i));
    YY_B_SP(i)=I3(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=I1(x1(i),y1(i));
    YY_R_CZ(i)=I1(x1(i),y1(i)+1);
    XX_G_CZ(i)=I2(x1(i),y1(i));
    YY_G_CZ(i)=I2(x1(i),y1(i)+1);
    XX_B_CZ(i)=I3(x1(i),y1(i));
    YY_B_CZ(i)=I3(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=I1(x1(i),y1(i));
    YY_R_DJX(i)=I1(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=I2(x1(i),y1(i));
    YY_G_DJX(i)=I2(x1(i)+1,y1(i)+1);
    XX_B_DJX(i)=I3(x1(i),y1(i));
    YY_B_DJX(i)=I3(x1(i)+1,y1(i)+1);
end
%水平
R=figure('color',[1 1 1]);
scatter(XX_R_SP,YY_R_SP,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x+1,y)');
% title('原始图像R通道水平相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_original_Horizontal=frame2im(frame);
imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/lena/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/辣椒/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/狒狒/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/飞机/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/汽车/加密前R通道水平相关性点图.png','png');
% imwrite(R_original_Horizontal,'../相邻像素相关性分析图片/城堡/加密前R通道水平相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_SP,YY_G_SP,18,'filled');
% xlabel('G通道随机点像素灰度值');
% ylabel('与该点相邻水平方向像素灰度值');
% title('原始图像G通道水平相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_original_Horizontal=frame2im(frame);
% imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/lena/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/辣椒/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/狒狒/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/飞机/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/汽车/加密前G通道水平相关性点图.png','png');
% % imwrite(G_original_Horizontal,'../相邻像素相关性分析图片/城堡/加密前G通道水平相关性点图.png','png');
% 
% R=figure('color',[1 1 1]);
% scatter(XX_B_SP,YY_B_SP,18,'filled');
% xlabel('B通道随机点像素灰度值');
% ylabel('与该点相邻水平方向像素灰度值');
% title('原始图像B通道水平相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_original_Horizontal=frame2im(frame);
% imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/lena/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/辣椒/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/狒狒/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/飞机/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/汽车/加密前B通道水平相关性点图.png','png');
% % imwrite(B_original_Horizontal,'../相邻像素相关性分析图片/城堡/加密前B通道水平相关性点图.png','png');
%垂直
R=figure('color',[1 1 1]);
scatter(XX_R_CZ,YY_R_CZ,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x,y+1)');
% title('原始图像R通道垂直相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_original_Vertical=frame2im(frame);
imwrite(R_original_Vertical,'../相邻像素相关性分析图片/lena/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/辣椒/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/狒狒/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/飞机/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/汽车/加密前R通道垂直相关性点图.png','png');
% imwrite(R_original_Vertical,'../相邻像素相关性分析图片/城堡/加密前R通道垂直相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_CZ,YY_G_CZ,18,'filled');
% xlabel('G通道随机点像素灰度值');
% ylabel('与该点相邻垂直方向像素灰度值');
% title('原始图像G通道垂直相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_original_Vertical=frame2im(frame);
% imwrite(G_original_Vertical,'../相邻像素相关性分析图片/lena/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/辣椒/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/狒狒/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/飞机/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/汽车/加密前G通道垂直相关性点图.png','png');
% % imwrite(G_original_Vertical,'../相邻像素相关性分析图片/城堡/加密前G通道垂直相关性点图.png','png');
% 
% R=figure('color',[1 1 1]);
% scatter(XX_B_CZ,YY_B_CZ,18,'filled');
% xlabel('B通道随机点像素灰度值');
% ylabel('与该点相邻垂直方向像素灰度值');
% title('原始图像B通道垂直相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_original_Vertical=frame2im(frame);
% imwrite(B_original_Vertical,'../相邻像素相关性分析图片/lena/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/辣椒/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/狒狒/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/飞机/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/汽车/加密前B通道垂直相关性点图.png','png');
% % imwrite(B_original_Vertical,'../相邻像素相关性分析图片/城堡/加密前B通道垂直相关性点图.png','png');

%对角线
R=figure('color',[1 1 1]);
scatter(XX_R_DJX,YY_R_DJX,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x+1,y+1)');
% title('原始图像R通道对角线相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_original_Diagonal=frame2im(frame);
imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/lena/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/辣椒/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/狒狒/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/飞机/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/汽车/加密前R通道对角线相关性点图.png','png');
% imwrite(R_original_Diagonal,'../相邻像素相关性分析图片/城堡/加密前R通道对角线相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_DJX,YY_G_DJX,18,'filled');
% xlabel('G通道随机点像素灰度值');
% ylabel('与该点相邻对角线方向像素灰度值');
% title('原始图像G通道对角线相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_original_Diagonal=frame2im(frame);
% imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/lena/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/辣椒/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/狒狒/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/飞机/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/汽车/加密前G通道对角线相关性点图.png','png');
% % imwrite(G_original_Diagonal,'../相邻像素相关性分析图片/城堡/加密前G通道对角线相关性点图.png','png');
% 
% R=figure('color',[1 1 1]);
% scatter(XX_B_DJX,YY_B_DJX,18,'filled');
% xlabel('B通道随机点像素灰度值');
% ylabel('与该点相邻对角线方向像素灰度值');
% title('原始图像B通道对角线相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_original_Diagonal=frame2im(frame);
% imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/lena/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/辣椒/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/狒狒/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/飞机/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/汽车/加密前B通道对角线相关性点图.png','png');
% % imwrite(B_original_Diagonal,'../相邻像素相关性分析图片/城堡/加密前B通道对角线相关性点图.png','png');

%R通道
EX1_R=0;EY1_SP_R=0;DX1_R=0;DY1_SP_R=0;COVXY1_SP_R=0;    %计算水平相关性时需要的变量
EY1_CZ_R=0;DY1_CZ_R=0;COVXY1_CZ_R=0;                %垂直
EY1_DJX_R=0;DY1_DJX_R=0;COVXY1_DJX_R=0;             %对角线
%G通道
EX1_G=0;EY1_SP_G=0;DX1_G=0;DY1_SP_G=0;COVXY1_SP_G=0;
EY1_CZ_G=0;DY1_CZ_G=0;COVXY1_CZ_G=0;
EY1_DJX_G=0;DY1_DJX_G=0;COVXY1_DJX_G=0;
%B通道
EX1_B=0;EY1_SP_B=0;DX1_B=0;DY1_SP_B=0;COVXY1_SP_B=0;
EY1_CZ_B=0;DY1_CZ_B=0;COVXY1_CZ_B=0;
EY1_DJX_B=0;DY1_DJX_B=0;COVXY1_DJX_B=0;

I1=double(I1);%将I1转换成双精度浮点类型
I2=double(I2);%将I2转换成双精度浮点类型
I3=double(I3);%将I3转换成双精度浮点类型
for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX1表示
    EX1_R=EX1_R+I1(x1(i),y1(i)); 
    EX1_G=EX1_G+I2(x1(i),y1(i)); 
    EX1_B=EX1_B+I3(x1(i),y1(i)); 
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY1_SP、EY1_CZ、EY1_DJX
    %R通道
    EY1_SP_R=EY1_SP_R+I1(x1(i),y1(i)+1);
    EY1_CZ_R=EY1_CZ_R+I1(x1(i)+1,y1(i));
    EY1_DJX_R=EY1_DJX_R+I1(x1(i)+1,y1(i)+1);
    %G通道
    EY1_SP_G=EY1_SP_G+I2(x1(i),y1(i)+1);
    EY1_CZ_G=EY1_CZ_G+I2(x1(i)+1,y1(i));
    EY1_DJX_G=EY1_DJX_G+I2(x1(i)+1,y1(i)+1);
    %B通道
    EY1_SP_B=EY1_SP_B+I3(x1(i),y1(i)+1);
    EY1_CZ_B=EY1_CZ_B+I3(x1(i)+1,y1(i));
    EY1_DJX_B=EY1_DJX_B+I3(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
% R通道
EX1_R=EX1_R/NN;
EY1_SP_R=EY1_SP_R/NN;
EY1_CZ_R=EY1_CZ_R/NN;
EY1_DJX_R=EY1_DJX_R/NN;
% G通道
EX1_G=EX1_G/NN;
EY1_SP_G=EY1_SP_G/NN;
EY1_CZ_G=EY1_CZ_G/NN;
EY1_DJX_G=EY1_DJX_G/NN;
% B通道
EX1_B=EX1_B/NN;
EY1_SP_B=EY1_SP_B/NN;
EY1_CZ_B=EY1_CZ_B/NN;
EY1_DJX_B=EY1_DJX_B/NN;
for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX1_R=DX1_R+(I1(x1(i),y1(i))-EX1_R)^2;
    DX1_G=DX1_G+(I2(x1(i),y1(i))-EX1_G)^2;
    DX1_B=DX1_B+(I3(x1(i),y1(i))-EX1_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    %R通道
    DY1_SP_R=DY1_SP_R+(I1(x1(i),y1(i)+1)-EY1_SP_R)^2;
    DY1_CZ_R=DY1_CZ_R+(I1(x1(i)+1,y1(i))-EY1_CZ_R)^2;
    DY1_DJX_R=DY1_DJX_R+(I1(x1(i)+1,y1(i)+1)-EY1_DJX_R)^2;
    %G通道
    DY1_SP_G=DY1_SP_G+(I2(x1(i),y1(i)+1)-EY1_SP_G)^2;
    DY1_CZ_G=DY1_CZ_G+(I2(x1(i)+1,y1(i))-EY1_CZ_G)^2;
    DY1_DJX_G=DY1_DJX_G+(I2(x1(i)+1,y1(i)+1)-EY1_DJX_G)^2;
    %B通道
    DY1_SP_B=DY1_SP_B+(I3(x1(i),y1(i)+1)-EY1_SP_B)^2;
    DY1_CZ_B=DY1_CZ_B+(I3(x1(i)+1,y1(i))-EY1_CZ_B)^2;
    DY1_DJX_B=DY1_DJX_B+(I3(x1(i)+1,y1(i)+1)-EY1_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %R通道
    COVXY1_SP_R=COVXY1_SP_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i),y1(i)+1)-EY1_SP_R);
    COVXY1_CZ_R=COVXY1_CZ_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i)+1,y1(i))-EY1_CZ_R);
    COVXY1_DJX_R=COVXY1_DJX_R+(I1(x1(i),y1(i))-EX1_R)*(I1(x1(i)+1,y1(i)+1)-EY1_DJX_R);
    %G通道
    COVXY1_SP_G=COVXY1_SP_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i),y1(i)+1)-EY1_SP_G);
    COVXY1_CZ_G=COVXY1_CZ_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i)+1,y1(i))-EY1_CZ_G);
    COVXY1_DJX_G=COVXY1_DJX_G+(I2(x1(i),y1(i))-EX1_G)*(I2(x1(i)+1,y1(i)+1)-EY1_DJX_G);
    %B通道
    COVXY1_SP_B=COVXY1_SP_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i),y1(i)+1)-EY1_SP_B);
    COVXY1_CZ_B=COVXY1_CZ_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i)+1,y1(i))-EY1_CZ_B);
    COVXY1_DJX_B=COVXY1_DJX_B+(I3(x1(i),y1(i))-EX1_B)*(I3(x1(i)+1,y1(i)+1)-EY1_DJX_B);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
DX1_R=DX1_R/NN;
DY1_SP_R=DY1_SP_R/NN;
DY1_CZ_R=DY1_CZ_R/NN;
DY1_DJX_R=DY1_DJX_R/NN;
COVXY1_SP_R=COVXY1_SP_R/NN;
COVXY1_CZ_R=COVXY1_CZ_R/NN;
COVXY1_DJX_R=COVXY1_DJX_R/NN;
%G通道
DX1_G=DX1_G/NN;
DY1_SP_G=DY1_SP_G/NN;
DY1_CZ_G=DY1_CZ_G/NN;
DY1_DJX_G=DY1_DJX_G/NN;
COVXY1_SP_G=COVXY1_SP_G/NN;
COVXY1_CZ_G=COVXY1_CZ_G/NN;
COVXY1_DJX_G=COVXY1_DJX_G/NN;
%B通道
DX1_B=DX1_B/NN;
DY1_SP_B=DY1_SP_B/NN;
DY1_CZ_B=DY1_CZ_B/NN;
DY1_DJX_B=DY1_DJX_B/NN;
COVXY1_SP_B=COVXY1_SP_B/NN;
COVXY1_CZ_B=COVXY1_CZ_B/NN;
COVXY1_DJX_B=COVXY1_DJX_B/NN;
%水平、垂直、对角线的相关性
%R通道
RXY1_SP_R=COVXY1_SP_R/sqrt(DX1_R*DY1_SP_R);
RXY1_CZ_R=COVXY1_CZ_R/sqrt(DX1_R*DY1_CZ_R);
RXY1_DJX_R=COVXY1_DJX_R/sqrt(DX1_R*DY1_DJX_R);
%G通道
RXY1_SP_G=COVXY1_SP_G/sqrt(DX1_G*DY1_SP_G);
RXY1_CZ_G=COVXY1_CZ_G/sqrt(DX1_G*DY1_CZ_G);
RXY1_DJX_G=COVXY1_DJX_G/sqrt(DX1_G*DY1_DJX_G);
%B通道
RXY1_SP_B=COVXY1_SP_B/sqrt(DX1_B*DY1_SP_B);
RXY1_CZ_B=COVXY1_CZ_B/sqrt(DX1_B*DY1_CZ_B);
RXY1_DJX_B=COVXY1_DJX_B/sqrt(DX1_B*DY1_DJX_B);

%% 加密图像相邻图像相关性分析
%{
先随机在0~M-1行和0~N-1列选中1000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
%相关性曲线
%水平
XX_R_SP=zeros(1,NN);YY_R_SP=zeros(1,NN);  %预分配内存
XX_G_SP=zeros(1,NN);YY_G_SP=zeros(1,NN);
XX_B_SP=zeros(1,NN);YY_B_SP=zeros(1,NN);
%垂直
XX_R_CZ=zeros(1,NN);YY_R_CZ=zeros(1,NN);  %预分配内存
XX_G_CZ=zeros(1,NN);YY_G_CZ=zeros(1,NN);
XX_B_CZ=zeros(1,NN);YY_B_CZ=zeros(1,NN);
%对角线
XX_R_DJX=zeros(1,NN);YY_R_DJX=zeros(1,NN);  %预分配内存
XX_G_DJX=zeros(1,NN);YY_G_DJX=zeros(1,NN);
XX_B_DJX=zeros(1,NN);YY_B_DJX=zeros(1,NN);
for i=1:NN
    %水平
    XX_R_SP(i)=Q_R(x1(i),y1(i));
    YY_R_SP(i)=Q_R(x1(i)+1,y1(i));
    XX_G_SP(i)=Q_G(x1(i),y1(i));
    YY_G_SP(i)=Q_G(x1(i)+1,y1(i));
    XX_B_SP(i)=Q_B(x1(i),y1(i));
    YY_B_SP(i)=Q_B(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=Q_R(x1(i),y1(i));
    YY_R_CZ(i)=Q_R(x1(i),y1(i)+1);
    XX_G_CZ(i)=Q_G(x1(i),y1(i));
    YY_G_CZ(i)=Q_G(x1(i),y1(i)+1);
    XX_B_CZ(i)=Q_B(x1(i),y1(i));
    YY_B_CZ(i)=Q_B(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=Q_R(x1(i),y1(i));
    YY_R_DJX(i)=Q_R(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=Q_G(x1(i),y1(i));
    YY_G_DJX(i)=Q_G(x1(i)+1,y1(i)+1);
    XX_B_DJX(i)=Q_B(x1(i),y1(i));
    YY_B_DJX(i)=Q_B(x1(i)+1,y1(i)+1);
end
%水平
R=figure('color',[1 1 1]);
scatter(XX_R_SP,YY_R_SP,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x+1,y)');
% title('加密图像R通道水平相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_aftet_Horizontal=frame2im(frame);
imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/lena/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/辣椒/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/狒狒/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/飞机/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/汽车/加密后R通道水平相关性点图.png','png');
% imwrite(R_aftet_Horizontal,'../相邻像素相关性分析图片/城堡/加密后R通道水平相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_SP,YY_G_SP,18,'filled');
% xlabel('G通道随机点像素灰度值');
% ylabel('与该点相邻水平方向像素灰度值');
% title('加密图像G通道水平相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_aftet_Horizontal=frame2im(frame);
% imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/lena/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/辣椒/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/狒狒/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/飞机/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/汽车/加密后G通道水平相关性点图.png','png');
% % imwrite(G_aftet_Horizontal,'../相邻像素相关性分析图片/城堡/加密后G通道水平相关性点图.png','png');
% 
% R=figure('color',[1 1 1]);
% scatter(XX_B_SP,YY_B_SP,18,'filled');
% xlabel('B通道随机点像素灰度值');
% ylabel('与该点相邻水平方向像素灰度值');
% title('加密图像B通道水平相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_aftet_Horizontal=frame2im(frame);
% imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/lena/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/辣椒/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/狒狒/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/飞机/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/4.1.05房子/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/4.2.06小船/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/汽车/加密后B通道水平相关性点图.png','png');
% % imwrite(B_aftet_Horizontal,'../相邻像素相关性分析图片/城堡/加密后B通道水平相关性点图.png','png');

%垂直
R=figure('color',[1 1 1]);
scatter(XX_R_CZ,YY_R_CZ,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x,y+1)');
% title('加密图像R通道垂直相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_aftet_Vertical=frame2im(frame);
imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/lena/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/辣椒/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/狒狒/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/飞机/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/汽车/加密后R通道垂直相关性点图.png','png');
% imwrite(R_aftet_Vertical,'../相邻像素相关性分析图片/城堡/加密后R通道垂直相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_CZ,YY_G_CZ,18,'filled');
% xlabel('G通道随机点像素灰度值');
% ylabel('与该点相邻垂直方向像素灰度值');
% title('加密图像G通道垂直相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_aftet_Vertical=frame2im(frame);
% imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/lena/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/辣椒/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/狒狒/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/飞机/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/汽车/加密后G通道垂直相关性点图.png','png');
% % imwrite(G_aftet_Vertical,'../相邻像素相关性分析图片/城堡/加密后G通道垂直相关性点图.png','png');
% 
% R=figure('color',[1 1 1]);
% scatter(XX_B_CZ,YY_B_CZ,18,'filled');
% xlabel('B通道随机点像素灰度值');
% ylabel('与该点相邻垂直方向像素灰度值');
% title('加密图像B通道垂直相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_aftet_Vertical=frame2im(frame);
% imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/lena/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/辣椒/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/狒狒/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/飞机/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/4.1.05房子/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/4.2.06小船/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/汽车/加密后B通道垂直相关性点图.png','png');
% % imwrite(B_aftet_Vertical,'../相邻像素相关性分析图片/城堡/加密后B通道垂直相关性点图.png','png');

%对角线
R=figure('color',[1 1 1]);
scatter(XX_R_DJX,YY_R_DJX,18,'filled');
xlabel('Pixel value on location (x,y)');
ylabel('Pixel value on location (x+1,y+1)');
% title('加密图像R通道对角线相邻元素相关性点图');
axis([0 255,0 255]);
set(gca,'XTick',0:15:255);
set(gca,'YTick',0:15:255);
frame = getframe(R);
R_aftet_Diagonal=frame2im(frame);
imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/lena/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/辣椒/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/狒狒/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/飞机/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/汽车/加密后R通道对角线相关性点图.png','png');
% imwrite(R_aftet_Diagonal,'../相邻像素相关性分析图片/城堡/加密后R通道对角线相关性点图.png','png');

% R=figure('color',[1 1 1]);
% scatter(XX_G_DJX,YY_G_DJX,18,'filled');
% xlabel('Gray value of random pixels in channel G');
% ylabel('Gray value of pixel in diagonal direction adjacent to the point');
% % title('加密图像G通道对角线相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% G_aftet_Diagonal=frame2im(frame);
% imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/lena/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/辣椒/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/狒狒/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/飞机/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/汽车/加密后G通道对角线相关性点图.png','png');
% % imwrite(G_aftet_Diagonal,'../相邻像素相关性分析图片/城堡/加密后G通道对角线相关性点图.png','png');
% 
% % figure;scatter(XX_B_DJX,YY_B_DJX,18,'filled');xlabel('B通道随机点像素灰度值');ylabel('与该点相邻对角线方向像素灰度值');title('加密图像B通道对角线相邻元素相关性点图');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
% R=figure('color',[1 1 1]);
% scatter(XX_B_DJX,YY_B_DJX,18,'filled');
% xlabel('Gray value of random pixels in channel B');
% ylabel('Gray value of pixel in diagonal direction adjacent to the point');
% % title('加密图像B通道对角线相邻元素相关性点图');
% axis([0 255,0 255]);
% set(gca,'XTick',0:15:255);
% set(gca,'YTick',0:15:255);
% frame = getframe(R);
% B_aftet_Diagonal=frame2im(frame);
% imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/lena/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/辣椒/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/狒狒/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/飞机/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/4.1.05房子/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/4.2.06小船/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/汽车/加密后B通道对角线相关性点图.png','png');
% % imwrite(B_aftet_Diagonal,'../相邻像素相关性分析图片/城堡/加密后B通道对角线相关性点图.png','png');

%R通道
Q_R=double(Q_R);
EX2_R=0;EY2_SP_R=0;DX2_R=0;DY2_SP_R=0;COVXY2_SP_R=0;    %水平
EY2_CZ_R=0;DY2_CZ_R=0;COVXY2_CZ_R=0;    %垂直
EY2_DJX_R=0;DY2_DJX_R=0;COVXY2_DJX_R=0;   %对角线
%G通道
Q_G=double(Q_G);
EX2_G=0;EY2_SP_G=0;DX2_G=0;DY2_SP_G=0;COVXY2_SP_G=0;    %水平
EY2_CZ_G=0;DY2_CZ_G=0;COVXY2_CZ_G=0;    %垂直
EY2_DJX_G=0;DY2_DJX_G=0;COVXY2_DJX_G=0;   %对角线
%B通道
Q_B=double(Q_B);
EX2_B=0;EY2_SP_B=0;DX2_B=0;DY2_SP_B=0;COVXY2_SP_B=0;    %水平
EY2_CZ_B=0;DY2_CZ_B=0;COVXY2_CZ_B=0;    %垂直
EY2_DJX_B=0;DY2_DJX_B=0;COVXY2_DJX_B=0;   %对角线
for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX2表示
    EX2_R=EX2_R+Q_R(x1(i),y1(i));
    EX2_G=EX2_G+Q_G(x1(i),y1(i));
    EX2_B=EX2_B+Q_B(x1(i),y1(i));
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY2_SP、EY2_CZ、EY2_DJX
    %R通道
    EY2_SP_R=EY2_SP_R+Q_R(x1(i),y1(i)+1);
    EY2_CZ_R=EY2_CZ_R+Q_R(x1(i)+1,y1(i));
    EY2_DJX_R=EY2_DJX_R+Q_R(x1(i)+1,y1(i)+1);
    %G通道
    EY2_SP_G=EY2_SP_G+Q_G(x1(i),y1(i)+1);
    EY2_CZ_G=EY2_CZ_G+Q_G(x1(i)+1,y1(i));
    EY2_DJX_G=EY2_DJX_G+Q_G(x1(i)+1,y1(i)+1);
    %B通道
    EY2_SP_B=EY2_SP_B+Q_B(x1(i),y1(i)+1);
    EY2_CZ_B=EY2_CZ_B+Q_B(x1(i)+1,y1(i));
    EY2_DJX_B=EY2_DJX_B+Q_B(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
EX2_R=EX2_R/NN;
EY2_SP_R=EY2_SP_R/NN;
EY2_CZ_R=EY2_CZ_R/NN;
EY2_DJX_R=EY2_DJX_R/NN;
%G通道
EX2_G=EX2_G/NN;
EY2_SP_G=EY2_SP_G/NN;
EY2_CZ_G=EY2_CZ_G/NN;
EY2_DJX_G=EY2_DJX_G/NN;
%B通道
EX2_B=EX2_B/NN;
EY2_SP_B=EY2_SP_B/NN;
EY2_CZ_B=EY2_CZ_B/NN;
EY2_DJX_B=EY2_DJX_B/NN;

for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX2表示
    DX2_R=DX2_R+(Q_R(x1(i),y1(i))-EX2_R)^2;
    DX2_G=DX2_G+(Q_G(x1(i),y1(i))-EX2_G)^2;
    DX2_B=DX2_B+(Q_B(x1(i),y1(i))-EX2_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY2_SP、DY2_CZ、DY2_DJX
    %R通道
    DY2_SP_R=DY2_SP_R+(Q_R(x1(i),y1(i)+1)-EY2_SP_R)^2;
    DY2_CZ_R=DY2_CZ_R+(Q_R(x1(i)+1,y1(i))-EY2_CZ_R)^2;
    DY2_DJX_R=DY2_DJX_R+(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R)^2;
    %G通道
    DY2_SP_G=DY2_SP_G+(Q_G(x1(i),y1(i)+1)-EY2_SP_G)^2;
    DY2_CZ_G=DY2_CZ_G+(Q_G(x1(i)+1,y1(i))-EY2_CZ_G)^2;
    DY2_DJX_G=DY2_DJX_G+(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G)^2;
    %B通道
    DY2_SP_B=DY2_SP_B+(Q_B(x1(i),y1(i)+1)-EY2_SP_B)^2;
    DY2_CZ_B=DY2_CZ_B+(Q_B(x1(i)+1,y1(i))-EY2_CZ_B)^2;
    DY2_DJX_B=DY2_DJX_B+(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %R通道
    COVXY2_SP_R=COVXY2_SP_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i),y1(i)+1)-EY2_SP_R);
    COVXY2_CZ_R=COVXY2_CZ_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i))-EY2_CZ_R);
    COVXY2_DJX_R=COVXY2_DJX_R+(Q_R(x1(i),y1(i))-EX2_R)*(Q_R(x1(i)+1,y1(i)+1)-EY2_DJX_R);
    %G通道
    COVXY2_SP_G=COVXY2_SP_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i),y1(i)+1)-EY2_SP_G);
    COVXY2_CZ_G=COVXY2_CZ_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i))-EY2_CZ_G);
    COVXY2_DJX_G=COVXY2_DJX_G+(Q_G(x1(i),y1(i))-EX2_G)*(Q_G(x1(i)+1,y1(i)+1)-EY2_DJX_G);
    %B通道
    COVXY2_SP_B=COVXY2_SP_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i),y1(i)+1)-EY2_SP_B);
    COVXY2_CZ_B=COVXY2_CZ_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i))-EY2_CZ_B);
    COVXY2_DJX_B=COVXY2_DJX_B+(Q_B(x1(i),y1(i))-EX2_B)*(Q_B(x1(i)+1,y1(i)+1)-EY2_DJX_B);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%R通道
DX2_R=DX2_R/NN;
DY2_SP_R=DY2_SP_R/NN;
DY2_CZ_R=DY2_CZ_R/NN;
DY2_DJX_R=DY2_DJX_R/NN;
COVXY2_SP_R=COVXY2_SP_R/NN;
COVXY2_CZ_R=COVXY2_CZ_R/NN;
COVXY2_DJX_R=COVXY2_DJX_R/NN;
%G通道
DX2_G=DX2_G/NN;
DY2_SP_G=DY2_SP_G/NN;
DY2_CZ_G=DY2_CZ_G/NN;
DY2_DJX_G=DY2_DJX_G/NN;
COVXY2_SP_G=COVXY2_SP_G/NN;
COVXY2_CZ_G=COVXY2_CZ_G/NN;
COVXY2_DJX_G=COVXY2_DJX_G/NN;
%B通道
DX2_B=DX2_B/NN;
DY2_SP_B=DY2_SP_B/NN;
DY2_CZ_B=DY2_CZ_B/NN;
DY2_DJX_B=DY2_DJX_B/NN;
COVXY2_SP_B=COVXY2_SP_B/NN;
COVXY2_CZ_B=COVXY2_CZ_B/NN;
COVXY2_DJX_B=COVXY2_DJX_B/NN;
%水平、垂直、对角线的相关性
%R通道
RXY2_SP_R=COVXY2_SP_R/sqrt(DX2_R*DY2_SP_R);
RXY2_CZ_R=COVXY2_CZ_R/sqrt(DX2_R*DY2_CZ_R);
RXY2_DJX_R=COVXY2_DJX_R/sqrt(DX2_R*DY2_DJX_R);
%G通道
RXY2_SP_G=COVXY2_SP_G/sqrt(DX2_G*DY2_SP_G);
RXY2_CZ_G=COVXY2_CZ_G/sqrt(DX2_G*DY2_CZ_G);
RXY2_DJX_G=COVXY2_DJX_G/sqrt(DX2_G*DY2_DJX_G);
%B通道
RXY2_SP_B=COVXY2_SP_B/sqrt(DX2_B*DY2_SP_B);
RXY2_CZ_B=COVXY2_CZ_B/sqrt(DX2_B*DY2_CZ_B);
RXY2_DJX_B=COVXY2_DJX_B/sqrt(DX2_B*DY2_DJX_B);

disp('R通道相关性：');
disp(['原始图片R通道相关性：','  水平相关性=',num2str(RXY1_SP_R),'    垂直相关性=',num2str(RXY1_CZ_R),'  对角线相关性=',num2str(RXY1_DJX_R)]);
disp(['加密图片R通道相关性：','  水平相关性=',num2str(RXY2_SP_R),'  垂直相关性=',num2str(RXY2_CZ_R),'  对角线相关性=',num2str(RXY2_DJX_R)]);
disp('G通道相关性：');
disp(['原始图片G通道相关性：','  水平相关性=',num2str(RXY1_SP_G),'   垂直相关性=',num2str(RXY1_CZ_G),'  对角线相关性=',num2str(RXY1_DJX_G)]);
disp(['加密图片G通道相关性：','  水平相关性=',num2str(RXY2_SP_G),' 垂直相关性=',num2str(RXY2_CZ_G),'  对角线相关性=',num2str(RXY2_DJX_G)]);
disp('B通道相关性：');
disp(['原始图片B通道相关性：','  水平相关性=',num2str(RXY1_SP_B),'   垂直相关性=',num2str(RXY1_CZ_B),'  对角线相关性=',num2str(RXY1_DJX_B)]);
disp(['加密图片B通道相关性：','  水平相关性=',num2str(RXY2_SP_B),'  垂直相关性=',num2str(RXY2_CZ_B),'  对角线相关性=',num2str(RXY2_DJX_B)]);