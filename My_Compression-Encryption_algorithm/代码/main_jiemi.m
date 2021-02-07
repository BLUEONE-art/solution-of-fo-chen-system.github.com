%% 联合彩色数字图像加密-压缩系统（解密）
%   @author:董昊
%   @date:2018.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../原始、加密、解密图片/lena/加密后的lena.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/lena/原始lena/加密后的lena.png','png');             %读取图像信息

% I=imread('../原始、加密、解密图片/辣椒/加密后的peppers.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/狒狒/加密后的baboon.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/飞机/加密后的airplane.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/4.1.05房子/加密后的4.1.05.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/4.2.06小船/加密后的4.2.06.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/汽车/加密后的house.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/城堡/加密后的kodim08.png','png');             %读取图像信息

% I=imread('../PSNR分析/lena/加密后CR=1的lena.png','png'); 
% I=imread('../PSNR分析/lena/加密后CR=0.875的lena.png','png'); 
% I=imread('../PSNR分析/lena/加密后CR=0.75的lena.png','png'); 
% I=imread('../PSNR分析/lena/加密后CR=0.5的lena.png','png'); 
load('Encryption_Key.mat');

U1=I(:,:,1);  
U2=I(:,:,2);  
U3=I(:,:,3);  

[M,N]=size(U1);                      %将图像的行列赋值给M,N

%密钥敏感度分析
t=KEY(1,1);    %分块大小

SUM=M*N;

%密钥敏感度分析
% u=3.9999000000000001;
u=KEY(1,2);
% u=3.78774868768;

xx0=KEY(1,10);
% xx0=0.3884000000000001;

xx1=KEY(1,11);
% xx1=0.4133000000000001;

ppx=zeros(1,M+1000);        %预分配内存
ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;
for i=1:M+999                 %进行M+999次循环，共得到M+1000点（包括初值）
    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %进行M+999次循环，共得到M+1000点（包括初值）
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %去除前1000点，获得更好的随机性
ppy=ppy(1001:length(ppy));

[~,Ux]=sort(ppx,'descend');
[~,Uy]=sort(ppy,'descend');

for i=N:-1:1
    temp = U1(:,i);
    U1(:,i) = U1(:,Uy(i));
    U1(:,Uy(i)) = temp;
    temp = U2(:,i);
    U2(:,i) = U2(:,Uy(i));
    U2(:,Uy(i)) = temp;
    temp = U3(:,i);
    U3(:,i) = U3(:,Uy(i));
    U3(:,Uy(i)) = temp;
end
for i=M:-1:1
    temp = U1(i,:);
    U1(i,:) = U1(Ux(i),:);
    U1(Ux(i),:) = temp;
    temp = U2(i,:);
    U2(i,:) = U2(Ux(i),:);
    U2(Ux(i),:) = temp;
    temp = U3(i,:);
    U3(i,:) = U3(Ux(i),:);
    U3(Ux(i),:) = temp;
end
% %% 2.产生Logistic混沌序列
% % u=3.990000000000001; %密钥敏感性测试  10^-15
% % u=3.99;%密钥：Logistic参数μ
% % x0=0.7067000000000001; %密钥敏感性测试  10^-16
% % x0=0.5475; %密钥：Logistic初值x0
%  
% x0=KEY(1,3);
% % x0=0.547623608200001;
% % x0=0.8477234627836;
% 
% p=zeros(1,SUM+1000);
% p(1)=x0;
% for i=1:SUM+999                        %进行SUM+999次循环，产生SUM+1000个数据
%     p(i+1)=u*p(i)*(1-p(i));
% end
% p=p(1001:length(p));
% 
% %% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R
% p=mod(round(p*10^4),256);
% R=reshape(p,N,M)';  %转成M行N列
% 
% %% 4.求解混沌方程
% % 求四个初值X0,Y0,Z0,H0
% r=(M/t)*(N/t);
% X0=KEY(1,4);
% % X0=0.49790000000001;
% % X0=0.3743264725;
% 
% Y0=KEY(1,5);
% % Y0=0.4276000000000001;
% % Y0=0.983982478264;
% 
% Z0=KEY(1,6);
% % Z0=0.70860000000001;
% % Z0=0.8324737824678;
% 
% H0=KEY(1,7);
% % H0=0.78070000000001;
% % H0=0.394738947;
% 
% %%超混沌Chen系统
% % A=chen_output(X0,Y0,Z0,H0,r);
% % X=A(:,1);
% % X=X(3002:length(X));
% % Y=A(:,2);
% % Y=Y(3002:length(Y));
% % Z=A(:,3);
% % Z=Z(3002:length(Z));
% % H=A(:,4);
% % H=H(3002:length(H));
% 
% %%分数阶超混沌Chen系统
% % h=0.01;%设置步长
% % NN=25000;%运行混沌系统得到的序列个数
% % q=0.95;%分数阶阶数
% % z0=[X0 Y0 Z0 H0];
% % [~,y]=FrataSim(h,NN,z0,q);
% % X=y(1,:);
% % X=X(5000:5000+r);        %去除前3001项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
% % Y=y(2,:);
% % Y=Y(5000:5000+r);
% % Z=y(3,:);
% % Z=Z(5000:5000+r);
% % H=y(4,:);
% % % H=H(5000:length(H));
% % H=H(5000:5000+r);
% 
% %针对超过512*512的更大幅图像，需要更多的混沌序列
% h=0.01;%设置步长
% NN=35000;%运行混沌系统得到的序列个数
% q=0.95;%分数阶阶数
% 
% z0=[X0 Y0 Z0 H0];
% [~,y]=FrataSim(h,NN,z0,q);
% X=y(1,:);
% X=X(5000:5000+r);        %去除前3001项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
% Y=y(2,:);
% Y=Y(5000:5000+r);
% Z=y(3,:);
% Z=Z(5000:5000+r);
% H=y(4,:);
% % H=H(5000:length(H));
% H=H(5000:5000+r);
% %% 5.DNA编码
% %X,Y分别决定I和R的DNA编码方式，有8种，1~8
% X=mod(round(X*10^4),8)+1;
% Y=mod(round(Y*10^4),8)+1;
% Z=mod(round(Z*10^4),4);
% Z(Z==0)=4;      %加减法互换
% Z(Z==1)=0;
% Z(Z==4)=1;
% H=mod(round(H*10^4),8)+1;
% e=N/t;
% for i=r:-1:2
%     Q1_Y=DNA_bian(fenkuai(t,U1,i),H(i));
%     Q1_Cb=DNA_bian(fenkuai(t,U2,i),H(i));
%     Q1_Cr=DNA_bian(fenkuai(t,U3,i),H(i));
%     
%     Q1_last_Y=DNA_bian(fenkuai(t,U1,i-1),H(i-1));
%     Q1_last_Cb=DNA_bian(fenkuai(t,U2,i-1),H(i-1));
%     Q1_last_Cr=DNA_bian(fenkuai(t,U3,i-1),H(i-1));
%     
%     Q2_Y=DNA_yunsuan(Q1_Y,Q1_last_Y,Z(i));        %扩散前
%     Q2_Cb=DNA_yunsuan(Q1_Cb,Q1_last_Cb,Z(i));
%     Q2_Cr=DNA_yunsuan(Q1_Cr,Q1_last_Cr,Z(i));
% 
%     Q3=DNA_bian(fenkuai(t,R,i),Y(i));
%     
%     Q4_Y=DNA_yunsuan(Q2_Y,Q3,Z(i));
%     Q4_Cb=DNA_yunsuan(Q2_Cb,Q3,Z(i));
%     Q4_Cr=DNA_yunsuan(Q2_Cr,Q3,Z(i));
%     
%     xx=floor(i/e)+1;
%     yy=mod(i,e);
%     if yy==0
%         xx=xx-1;
%         yy=e;
%     end
%     U1((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Y,X(i));
%     U2((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cb,X(i));
%     U3((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cr,X(i));
% end
% Q5_Y=DNA_bian(fenkuai(t,U1,1),H(1));
% Q5_Cb=DNA_bian(fenkuai(t,U2,1),H(1));
% Q5_Cr=DNA_bian(fenkuai(t,U3,1),H(1));
% 
% Q6=DNA_bian(fenkuai(t,R,1),Y(1));
% 
% Q7_Y=DNA_yunsuan(Q5_Y,Q6,Z(1));
% Q7_Cb=DNA_yunsuan(Q5_Cb,Q6,Z(1));
% Q7_Cr=DNA_yunsuan(Q5_Cr,Q6,Z(1));
% 
% U1(1:t,1:t)=DNA_jie(Q7_Y,X(1));
% U2(1:t,1:t)=DNA_jie(Q7_Cb,X(1));
% U3(1:t,1:t)=DNA_jie(Q7_Cr,X(1));
% 
% U1=double(U1);
% U2=double(U2);
% U3=double(U3);
% 
% %% 逆数据处理
% rate=4;     %压缩矩阵1存在的行列数，可作为密钥
% max_Multiple_first=KEY(1,12);%解密密钥之一
% % max_Multiple_first=9999;%解密密钥之一
% 
% max_Multiple=KEY(1,13);%解密密钥之一
% % max_Multiple=666;
% 
% load('Mask.mat');
% load('Coefficient_Matrix.mat');
% T1=zeros(t,t);  % 压缩矩阵
% T1(1,1)=max_Multiple_first;
% for k=2:rate
%     T1(1,k)=max_Multiple;
% %     T1(1,k+1:rate)=1;
% end
% for i=2:rate
%     for j=1:rate%+1-i
%         T1(i,j)=max_Multiple;
% %         T1(i,j+1:rate)=1;
%     end
% end
% % T1=T1.*YS1;%没有压缩的情况
% % T1=T1.*YS2;
% % T1(4,3)=1;
% % T1(4,4)=1;
% 
% load('Coefficient_Matrix.mat');
% 
% U1=U1.*Index_Matrix1;
% U2=U2.*Index_Matrix2;
% U3=U3.*Index_Matrix3;
% 
% U1=U1+High_frequency1;
% U2=U2+High_frequency2;
% U3=U3+High_frequency3;
% 
% fun = @(block_struct) block_struct.data./T1;         %Matlab建议用blockproc函数代替blkproc函数
% J1 = blockproc(U1,[t t],fun);                        %经过数据处理之后的矩阵
% J2 = blockproc(U2,[t t],fun);                        %经过数据处理之后的矩阵
% J3 = blockproc(U3,[t t],fun);                        %经过数据处理之后的矩阵
% %% 逆DCT处理
% T=dctmtx(t);        %8阶DCT变换矩阵
% fun = @(block_struct) T'*block_struct.data*T;       %Matlab建议用blockproc函数代替blkproc函数
% K1 = blockproc(J1,[t t],fun);              %DCT逆变换，重构图像
% K2 = blockproc(J2,[t t],fun);              %DCT逆变换，重构图像
% K3 = blockproc(J3,[t t],fun);              %DCT逆变换，重构图像
% 
% I_jiemi(:,:,1)=K1;
% I_jiemi(:,:,2)=K2;
% I_jiemi(:,:,3)=K3;
% 
% P1=im2uint8(K1);
% P2=im2uint8(K2);
% P3=im2uint8(K3);
% V1(:,:,1)=P1;
% V1(:,:,2)=P2;
% V1(:,:,3)=P3;
% %% 6、去除加密时补的零
% M1=KEY(1,8);   %加密时补零的参数，M1=mod(M,t);作为密钥
% N1=KEY(1,9);   %加密时补零的参数，N1=mod(N,t);作为密钥
% if M1~=0
%     I_jiemi=I_jiemi(1:M-t+M1,:,:);
% end
% if N1~=0
%     I_jiemi=I_jiemi(:,1:N-t+N1,:);
% end
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,1));
% %title('解密图片R通道直方图');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/lena/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/辣椒/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/狒狒/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/飞机/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.1.05房子/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.2.06小船/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/汽车/解密图片R通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/城堡/解密图片R通道直方图.png','png');
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,2));
% % title('解密图片G通道直方图');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/lena/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/辣椒/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/狒狒/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/飞机/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.1.05房子/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.2.06小船/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/汽车/解密图片G通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/城堡/解密图片G通道直方图.png','png');
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,3));
% % title('解密图片B通道直方图');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/lena/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/辣椒/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/狒狒/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/飞机/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.1.05房子/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/4.2.06小船/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/汽车/解密图片B通道直方图.png','png');
% % imwrite(B_jiamizhifangtu,'../原始、加密、解密图片/城堡/解密图片B通道直方图.png','png');
% 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/解密后的lena.png','png'); 
% 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/u+10^-16.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/X(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/Y(0)+10^-16.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/Z(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/H(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/x0+10^-15.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/x1+10^-16.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/x2+10^-16.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/p+10^-16.png','png'); 
% 
% % imwrite(I_jiemi,'../原始、加密、解密图片/lena/all_dontknow.png','png');
% 
% % imwrite(I_jiemi,'../原始、加密、解密图片/辣椒/解密后的peppers.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/狒狒/解密后的baboon.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/飞机/解密后的airplane.png','png'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/4.1.05房子/解密后的4.1.05.tiff','tiff'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/4.2.06小船/解密后的4.2.06.tiff','tiff'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/汽车/解密后的house.tiff','tiff'); 
% % imwrite(I_jiemi,'../原始、加密、解密图片/城堡/解密后的kodim08.png','png'); 
% 
% disp('您输入的解密密钥为：');
% disp(['密钥1：t=',num2str(t),'     密钥2：u=',num2str(u),'    密钥3：x0=',num2str(x0),'    密钥4：x(0)=',num2str(X0),'   密钥5：y(0)=',num2str(Y0) '   密钥6：z(0)=',num2str(Z0)]);
% disp(['密钥7：h(0)=',num2str(H0),'   密钥8：M1=',num2str(M1),'   密钥9：N1=',num2str(N1),'   密钥10：xx0=',num2str(xx0),'   密钥11：xx1=',num2str(xx1) '   密钥12：max_Multiple_first=',num2str(max_Multiple_first) '   密钥13：max_Multiple=',num2str(max_Multiple)]);
% disp('解密完成'); 
% figure;imshow(I_jiemi);
% title('解密后图片');
% % figure;imshow(V1);
% % title('解密后图片2');