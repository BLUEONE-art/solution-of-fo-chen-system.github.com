%% 抗噪声攻击分析
%   @author:董昊
%   @date:2020.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../原始、加密、解密图片/lena/加密后的lena.png','png');           %读取图像信息
% I=imread('../原始、加密、解密图片/辣椒/加密后的peppers.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/狒狒/加密后的baboon.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/飞机/加密后的airplane.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/4.1.05房子/加密后的4.1.05.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/4.2.06小船/加密后的4.2.06.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/汽车/加密后的house.png','png');             %读取图像信息
% I=imread('../原始、加密、解密图片/城堡/加密后的kodim08.png','png');             %读取图像信息
load('Encryption_Key.mat');
load('Mask.mat');
U1=I(:,:,1);  
U2=I(:,:,2);  
U3=I(:,:,3);  

[M,N]=size(U1);                      %将图像的行列赋值给M,N
t=KEY(1,1);    %分块大小
SUM=M*N;
% u=3.9999;
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
I(:,:,1)=U1;
I(:,:,2)=U2;
I(:,:,3)=U3;
%% 2.产生Logistic混沌序列
% u=3.990000000000001; %密钥敏感性测试  10^-15
% u=3.99;%密钥：Logistic参数μ
% x0=0.7067000000000001; %密钥敏感性测试  10^-16
% x0=0.5475; %密钥：Logistic初值x0
 
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
%%超混沌Chen系统
% A=chen_output(X0,Y0,Z0,H0,r);
% X=A(:,1);
% X=X(3002:length(X));
% Y=A(:,2);
% Y=Y(3002:length(Y));
% Z=A(:,3);
% Z=Z(3002:length(Z));
% H=A(:,4);
% H=H(3002:length(H));

%%分数阶超混沌Chen系统
h=0.01;%设置步长
NN=25000;%运行混沌系统得到的序列个数
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
% H=H(5000:length(H));
H=H(5000:5000+r);

 %X,Y分别决定I和R的DNA编码方式，有8种，1~8
X=mod(round(X*10^4),8)+1;
Y=mod(round(Y*10^4),8)+1;
Z=mod(round(Z*10^4),4);
Z(Z==0)=4;
Z(Z==1)=0;
Z(Z==4)=1;
H=mod(round(H*10^4),8)+1;
e=N/t;
%% 图像质量评价
YY=imread('../原始、加密、解密图片/lena/lena.png','png');        %读取图像信息
YY=double(YY);
Y1=YY(:,:,1);        %R
Y2=YY(:,:,2);        %G
Y3=YY(:,:,3);        %B
MSE_R=zeros(1,21);MSE_G=zeros(1,21);MSE_B=zeros(1,21);
j=0;        %数组下标
for i=0:5:5
    P = imnoise(I, 'gaussian', 0, i^2/255^2);  %加入高斯白噪声
    imwrite(P,'../加入高斯噪声解密结果/lena/高斯噪声方差为5时的加密lena.png','png'); 
% for i=0:0.05:0.05
%     P = imnoise(I,'salt & pepper',i);
    U1=P(:,:,1);     %R通道
    U2=P(:,:,2);     %G通道
    U3=P(:,:,3);     %B通道
    j=j+1;      %数组下标

    %% DNA编码（解密）
    for ii=r:-1:2
        Q1_Y=DNA_bian(fenkuai(t,U1,ii),H(ii));
        Q1_Cb=DNA_bian(fenkuai(t,U2,ii),H(ii));
        Q1_Cr=DNA_bian(fenkuai(t,U3,ii),H(ii));

        Q1_last_Y=DNA_bian(fenkuai(t,U1,ii-1),H(ii-1));
        Q1_last_Cb=DNA_bian(fenkuai(t,U2,ii-1),H(ii-1));
        Q1_last_Cr=DNA_bian(fenkuai(t,U3,ii-1),H(ii-1));

        Q2_Y=DNA_yunsuan(Q1_Y,Q1_last_Y,Z(ii));        %扩散前
        Q2_Cb=DNA_yunsuan(Q1_Cb,Q1_last_Cb,Z(ii));
        Q2_Cr=DNA_yunsuan(Q1_Cr,Q1_last_Cr,Z(ii));

        Q3=DNA_bian(fenkuai(t,R,ii),Y(ii));

        Q4_Y=DNA_yunsuan(Q2_Y,Q3,Z(ii));
        Q4_Cb=DNA_yunsuan(Q2_Cb,Q3,Z(ii));
        Q4_Cr=DNA_yunsuan(Q2_Cr,Q3,Z(ii));

        xx=floor(ii/e)+1;
        yy=mod(ii,e);
        if yy==0
            xx=xx-1;
            yy=e;
        end
        U1((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Y,X(ii));
        U2((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cb,X(ii));
        U3((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_Cr,X(ii));
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
    rate=4;     %压缩矩阵1存在的行列数，可作为密钥
    max_Multiple_first=KEY(1,12);%解密密钥之一
    max_Multiple=KEY(1,13);%解密密钥之一
    load('Coefficient_Matrix.mat');
    T1=zeros(t,t);  % 压缩矩阵
    T1(1,1)=max_Multiple_first;
    for k=2:rate
        T1(1,k)=max_Multiple;
    %     T1(1,k+1:rate)=1;
    end
    for a=2:rate
        for b=1:rate
            T1(a,b)=max_Multiple;
%             T1(a,b+1:rate)=1;
        end
    end
%     T1=T1.*YS1;%没有压缩的情况
    % T1=T1.*YS2;
    % T1(4,3)=1;
    % T1(4,4)=1;

    load('Coefficient_Matrix.mat');
    U1=U1.*Index_Matrix1;
    U2=U2.*Index_Matrix2;
    U3=U3.*Index_Matrix3;
    U1=U1+High_frequency1;
    U2=U2+High_frequency2;
    U3=U3+High_frequency3;
    fun = @(block_struct) block_struct.data./T1;         %Matlab建议用blockproc函数代替blkproc函数
    J1 = blockproc(U1,[t t],fun);                        %经过数据处理之后的矩阵
    J2 = blockproc(U2,[t t],fun);                        %经过数据处理之后的矩阵
    J3 = blockproc(U3,[t t],fun);                        %经过数据处理之后的矩阵
    %% 逆DCT处理
    T=dctmtx(t);        %8阶DCT变换矩阵
    fun = @(block_struct) T'*block_struct.data*T;       %Matlab建议用blockproc函数代替blkproc函数
    K1 = blockproc(J1,[t t],fun);              %DCT逆变换，重构图像
    K2 = blockproc(J2,[t t],fun);              %DCT逆变换，重构图像
    K3 = blockproc(J3,[t t],fun);              %DCT逆变换，重构图像

%     I_jiemi(:,:,1)=K1;
%     I_jiemi(:,:,2)=K2;
%     I_jiemi(:,:,3)=K3;

    K1=im2uint8(K1);
    K2=im2uint8(K2);
    K3=im2uint8(K3);
    
    K1=double(K1);
    K2=double(K2);
    K3=double(K3);
    
    %% 6、去除加密时补的零
    M1=KEY(1,8);   %加密时补零的参数，M1=mod(M,t);作为密钥
    N1=KEY(1,9);   %加密时补零的参数，N1=mod(N,t);作为密钥
% %     if M1~=0
% %         I_jiemi=I_jiemi(1:M-t+M1,:,:);
% %     end
% %     if N1~=0
% %         I_jiemi=I_jiemi(:,1:N-t+N1,:);
% %     end
    %去除加密时补的零
    if M1~=0
        K1=K1(1:M-t+M1,:);
        K2=K2(1:M-t+M1,:);
        K3=K3(1:M-t+M1,:);
    end
    if N1~=0
        K1=K1(:,1:N-t+N1);
        K2=K2(:,1:N-t+N1);
        K3=K3(:,1:N-t+N1);
    end
    [MM,NN]=size(K1);     %重新获得解密后的图片大小
    for m=1:MM
        for n=1:NN
            MSE_R(j)=MSE_R(j)+(Y1(m,n)-K1(m,n))^2;       %R通道MSE
            MSE_G(j)=MSE_G(j)+(Y2(m,n)-K2(m,n))^2;       %G通道MSE
            MSE_B(j)=MSE_B(j)+(Y3(m,n)-K3(m,n))^2;       %B通道MSE
        end
    end

    I_jiemi(:,:,1)=uint8(K1);
    I_jiemi(:,:,2)=uint8(K2);
    I_jiemi(:,:,3)=uint8(K3);
%     I_jiemi(:,:,1)=K1;
%     I_jiemi(:,:,2)=K2;
%     I_jiemi(:,:,3)=K3;
% 
    figure;imshow(I_jiemi);title(['高斯噪声强度为',num2str(i),'时的解密图像']);
    imwrite(I_jiemi,'../加入高斯噪声解密结果/lena/高斯噪声方差为5时的lena解密.png','png'); 
%     figure;imshow(I_jiemi);title(['椒盐噪声密度为',num2str(i),'时的解密图像']);
%     imwrite(I_jiemi,'../加入椒盐噪声解密结果/lena/高斯噪声方差为5时的lena解密.png','png'); 
end
%噪声功率-MSE
MSE_R=MSE_R./SUM;
MSE_G=MSE_G./SUM;
MSE_B=MSE_B./SUM;
%峰值信噪比-PSNR
PSNR_R=10*log10((255^2)./MSE_R);
PSNR_G=10*log10((255^2)./MSE_G);
PSNR_B=10*log10((255^2)./MSE_B);
% 绘图，噪声功率-MSE、峰值信噪比-PSNR
% set(gca,'xtick', X);
% X=0:0.05:1;
% X=0:5:40;
% figure;plot(X,MSE_R);xlabel('椒盐噪声密度');ylabel('均方误差MSE');title('R通道：椒盐噪声密度-均方误差MSE曲线图');
% figure;plot(X,MSE_G);xlabel('椒盐噪声密度');ylabel('均方误差MSE');title('G通道：椒盐噪声密度-均方误差MSE曲线图');
% figure;plot(X,MSE_B);xlabel('椒盐噪声密度');ylabel('均方误差MSE');title('B通道：椒盐噪声密度-均方误差MSE曲线图');
% figure;plot(X,PSNR_R);xlabel('椒盐噪声密度');ylabel('峰值信噪比PSNR（dB）');title('R通道：椒盐噪声密度-峰值信噪比PSNR曲线图');
% figure;plot(X,PSNR_G);xlabel('椒盐噪声密度');ylabel('峰值信噪比PSNR（dB）');title('G通道：椒盐噪声密度-峰值信噪比PSNR曲线图');
% figure;plot(X,PSNR_B);xlabel('椒盐噪声密度');ylabel('峰值信噪比PSNR（dB）');title('B通道：椒盐噪声密度-峰值信噪比PSNR曲线图');