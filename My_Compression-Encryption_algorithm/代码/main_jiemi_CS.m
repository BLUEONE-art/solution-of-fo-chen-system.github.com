%% 联合彩色数字图像加密-压缩系统（解密）
%   @author:董昊
%   @date:2018.08.07
%-------------------------------------------------------------------------------------------------------%
clear;clc;
O=imread('../原始、加密、解密图片/lena_CS/lena.png','png');             %读取原始图像信息，用于计算PSNR
% O=double(O);
% I=imread('../原始、加密、解密图片/lena_CS/加密后的lena.png','png');             %读取图像信息

% I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.875的lena.png','png'); 
% I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.8125的lena.png','png'); 
% I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.75的lena.png','png'); 
% I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.5的lena.png','png'); 
% I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.4的lena.png','png'); 
I=imread('../原始、加密、解密图片/lena_CS/加密后CR=0.2的lena.png','png'); 


load('Encryption_Key.mat');

U1=I(:,:,1);  
U2=I(:,:,2);  
U3=I(:,:,3);  

K1 = O(:,:,1);
K2 = O(:,:,2);
K3 = O(:,:,3);

[M,N]=size(U1);                      %将图像的行列赋值给M,N
t=KEY(1,1);    %分块大小
SUM=M*N;
u=KEY(1,2);
xx0=KEY(1,10);
xx1=KEY(1,11);
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

%% 2.产生Logistic混沌序列
x0=KEY(1,3);

p=zeros(1,SUM+1000);
p(1)=x0;
for i=1:SUM+999                        %进行SUM+999次循环，产生SUM+1000个数据
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));

%% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R
p=mod(round(p*10^4),256);
R=reshape(p,N,M)';  %转成M行N列

%% 4.求解混沌方程
% 求四个初值X0,Y0,Z0,H0
r=(M/t)*(N/t);
X0=KEY(1,4);
Y0=KEY(1,5);
Z0=KEY(1,6);
H0=KEY(1,7);

%% 分数阶超混沌Chen系统
%针对超过512*512的更大幅图像，需要更多的混沌序列
h=0.01;%设置步长
NN=35000;%运行混沌系统得到的序列个数
q=0.95;%分数阶阶数

z0=[X0 Y0 Z0 H0];
[~,y]=FrataSim(h,NN,z0,q);
X=y(1,:);
X=X(5000:5000+r);        %去除前3001项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
Y=y(2,:);
Y=Y(5000:5000+r);
Z=y(3,:);
Z=Z(5000:5000+r);
H=y(4,:);
H=H(5000:5000+r);
%% 5.DNA编码
%X,Y分别决定I和R的DNA编码方式，有8种，1~8
X=mod(round(X*10^4),8)+1;
Y=mod(round(Y*10^4),8)+1;
Z=mod(round(Z*10^4),4);
Z(Z==0)=4;      %加减法互换
Z(Z==1)=0;
Z(Z==4)=1;
H=mod(round(H*10^4),8)+1;
e=N/t;
for i=r:-1:2
    Q1_Y=DNA_bian(fenkuai(t,U1,i),H(i));
    Q1_Cb=DNA_bian(fenkuai(t,U2,i),H(i));
    Q1_Cr=DNA_bian(fenkuai(t,U3,i),H(i));
    
    Q1_last_Y=DNA_bian(fenkuai(t,U1,i-1),H(i-1));
    Q1_last_Cb=DNA_bian(fenkuai(t,U2,i-1),H(i-1));
    Q1_last_Cr=DNA_bian(fenkuai(t,U3,i-1),H(i-1));
    
    Q2_Y=DNA_yunsuan(Q1_Y,Q1_last_Y,Z(i));        %扩散前
    Q2_Cb=DNA_yunsuan(Q1_Cb,Q1_last_Cb,Z(i));
    Q2_Cr=DNA_yunsuan(Q1_Cr,Q1_last_Cr,Z(i));

    Q3=DNA_bian(fenkuai(t,R,i),Y(i));
    
    Q4_Y=DNA_yunsuan(Q2_Y,Q3,Z(i));
    Q4_Cb=DNA_yunsuan(Q2_Cb,Q3,Z(i));
    Q4_Cr=DNA_yunsuan(Q2_Cr,Q3,Z(i));
    
    xx=floor(i/e)+1;
    yy=mod(i,e);
    if yy==0
        xx=xx-1;
        yy=e;
    end
    U1((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Y,X(i));
    U2((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cb,X(i));
    U3((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cr,X(i));
end
Q5_Y=DNA_bian(fenkuai(t,U1,1),H(1));
Q5_Cb=DNA_bian(fenkuai(t,U2,1),H(1));
Q5_Cr=DNA_bian(fenkuai(t,U3,1),H(1));

Q6=DNA_bian(fenkuai(t,R,1),Y(1));

Q7_Y=DNA_yunsuan(Q5_Y,Q6,Z(1));
Q7_Cb=DNA_yunsuan(Q5_Cb,Q6,Z(1));
Q7_Cr=DNA_yunsuan(Q5_Cr,Q6,Z(1));

U1(1:t,1:t)=DNA_jie(Q7_Y,X(1));
U2(1:t,1:t)=DNA_jie(Q7_Cb,X(1));
U3(1:t,1:t)=DNA_jie(Q7_Cr,X(1));

U1=double(U1);
U2=double(U2);
U3=double(U3);

%% 逆数据处理
load('Coefficient_Matrix.mat');

U1=U1.*Index_Matrix1;
U2=U2.*Index_Matrix2;
U3=U3.*Index_Matrix3;

U1=U1+High_frequency1;
U2=U2+High_frequency2;
U3=U3+High_frequency3;

J1 = U1;                        %经过数据处理之后的矩阵
J2 = U2;                        %经过数据处理之后的矩阵
J3 = U3;                        %经过数据处理之后的矩阵

%%  OMP算法
%  小波变换矩阵生成
load('PartHadamardMatrix.mat'); 
load('Measure_MTX.mat');

% M2=256;      % CR = 100%
% M2=224;      % CR = 87.5%
% M2=208;      % CR = 81.25%
% M2=192;      % CR = 75%
% M2=128;      % CR = 50%
% M2=104;      % CR = 40%
% M2=52;       % CR = 20%

D1=zeros(256,N);  %  恢复矩阵
for i=1:N  %  列循环       
    rec=omp(J1(:,i),K,256);
    D1(:,i)=rec;
end
D2=zeros(256,N);  %  恢复矩阵
for i=1:N  %  列循环       
    rec=omp(J2(:,i),K,256);
    D2(:,i)=rec;
end
D3=zeros(256,N);  %  恢复矩阵
for i=1:N  %  列循环       
    rec=omp(J3(:,i),K,256);
    D3(:,i)=rec;
end

E1=ww'*sparse(D1)*ww;  %  小波反变换
E11=full(E1);
E111=255.*E11;

E2=ww'*sparse(D2)*ww;  %  小波反变换
E22=full(E2);
E222=255.*E22;

E3=ww'*sparse(D3)*ww;  %  小波反变换
E33=full(E3);
E333=255.*E33;

I_jiemi(:,:,1)=uint8(E111);
I_jiemi(:,:,2)=uint8(E222);
I_jiemi(:,:,3)=uint8(E333);

%% 6、去除加密时补的零
M1=KEY(1,8);   %加密时补零的参数，M1=mod(M,t);作为密钥
N1=KEY(1,9);   %加密时补零的参数，N1=mod(N,t);作为密钥
if M1~=0
    I_jiemi=I_jiemi(1:M-t+M1,:,:);
end
if N1~=0
    I_jiemi=I_jiemi(:,1:N-t+N1,:);
end

figure('color',[1 1 1]);
imhist(I_jiemi(:,:,1));
%title('解密图片R通道直方图');

figure('color',[1 1 1]);
imhist(I_jiemi(:,:,2));
% title('解密图片G通道直方图');

figure('color',[1 1 1]);
imhist(I_jiemi(:,:,3));
% title('解密图片B通道直方图');

% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后的lena.png','png');

% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.875的lena.png','png'); 
% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.8125的lena.png','png'); 
% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.75的lena.png','png'); 
% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.5的lena.png','png'); 
% imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.4的lena.png','png'); 
imwrite(I_jiemi,'../原始、加密、解密图片/lena_CS/解密后CR=0.2的lena.png','png'); 

%  误差(R_PSNR)
errorx1=sum(sum(abs(E111-double(K1)).^2));       %  MSE误差
psnr1 = 10*log10(255*255/(errorx1/M/N))   %  PSNR
%dNMSE1= nmse(uint8(K1),uint8(E11))   %NMSE误差

%  误差(G_PSNR)
errorx2=sum(sum(abs(E222-double(K2)).^2));       %  MSE误差
psnr2 = 10*log10(255*255/(errorx2/M/N))   %  PSNR
%dNMSE2= nmse(uint8(K2),uint8(E22))   %NMSE误差

%  误差(B_PSNR)
errorx3=sum(sum(abs(E333-double(K3)).^2));       %  MSE误差
psnr3 = 10*log10(255*255/(errorx3/M/N))   %  PSNR
%dNMSE3= nmse(uint8(K3),uint8(E33))   %NMSE误差

disp('您输入的解密密钥为：');
disp(['密钥1：t=',num2str(t),'     密钥2：u=',num2str(u),'    密钥3：x0=',num2str(x0),'    密钥4：x(0)=',num2str(X0),'   密钥5：y(0)=',num2str(Y0) '   密钥6：z(0)=',num2str(Z0)]);
disp(['密钥7：h(0)=',num2str(H0),'   密钥8：M1=',num2str(M1),'   密钥9：N1=',num2str(N1),'   密钥10：xx0=',num2str(xx0),'   密钥11：xx1=',num2str(xx1)]);
disp('解密完成'); 
disp(['R通道的PSNR=',num2str(psnr1)]); 
disp(['G通道的PSNR=',num2str(psnr2)]); 
disp(['B通道的PSNR=',num2str(psnr3)]); 
figure;imshow(I_jiemi);
title('解密后图片');