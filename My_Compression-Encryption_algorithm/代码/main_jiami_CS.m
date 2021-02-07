%% 联合彩色数字图像加密-压缩系统（基于CS加密）
%   @author:董昊
%   @date:2020.08.07
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../原始、加密、解密图片/lena/lena.png','png');         %读取图像信息，lena
figure;imshow(I);title('原始图片');
load('Mask.mat');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B

KEY=zeros(1,11);

Q=figure('color',[1 1 1]);
imhist(I1);%title('原始图片R通道直方图');
frame = getframe(Q);
R_zhifangtu=frame2im(frame);

W=figure('color',[1 1 1]);
imhist(I2);%title('原始图片G通道直方图');
frame = getframe(W);
G_zhifangtu=frame2im(frame);

E=figure('color',[1 1 1]);
imhist(I3);%title('原始图片B通道直方图');
frame = getframe(E);
B_zhifangtu=frame2im(frame);

[M,N]=size(I1);                      %将图像的行列赋值给M,N
t=4;    %分块大小,密钥1
%% 1.补零
%将图像的行列数都补成可以被t整除的数，t为分块的大小。
M1=mod(M,t);    %可作为固定密钥，以便解码时可以去除补上的0，密钥2
N1=mod(N,t);    %可作为固定密钥，以便解码时可以去除补上的0，密钥3
                %比如M=5，则M1为1，进行下面补零操作
if M1~=0
    I1(M+1:M+t-M1,:)=0;
    I2(M+1:M+t-M1,:)=0;
    I3(M+1:M+t-M1,:)=0;
end
if N1~=0
    I1(:,N+1:N+t-N1)=0;
    I2(:,N+1:N+t-N1)=0;
    I3(:,N+1:N+t-N1)=0;
end
[M,N]=size(I1);  %补零后的行数和列数


% SUM=M*N; %CR = 100%
% SUM=224*N; %CR = 87.5%
% SUM=M*208; %CR = 81.25%
% SUM=M*192; %CR = 75%
% SUM=M*128; %CR = 50%
% SUM=M*104; %CR = 40%
SUM=M*72; %CR = 20%


%% 2.产生Logistic混沌序列
u=3.9999;     %Logistic参数μ，自定为3.9999，作为密钥4
x0=(sum(I1(:))+sum(I2(:)))/(255*256*256*2);     %计算得出Logistic初值x0，作为密钥5，这里是根据原始RGB彩色空间图像进行算的密钥，只跟明文相关
% x0=floor(x0*10^10)/10^10;     %保留10位小数
x0=roundn(x0,-10);


% p=zeros(1,SUM+1000);            %预分配内存 CR = 100%
% p=zeros(1,SUM+1000);            %预分配内存 CR = 87.5%
% p=zeros(1,SUM+1000);            %预分配内存 CR = 81.25%
% p=zeros(1,SUM+1000);            %预分配内存 CR = 75%
% p=zeros(1,SUM+1000);            %预分配内存 CR = 50%
% p=zeros(1,SUM+1000);            %预分配内存 CR = 40%
p=zeros(1,SUM+1000);            %预分配内存 CR = 20%


p(1)=x0;
for i=1:SUM+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));            %去除前1000点，获得更好的随机性

%% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R,这里要根据压缩率进行修改！
p=mod(round(p*10^4),256); 


% R=reshape(p,N,M)';  %转成M行N列的随机矩阵R CR = 100%
% R=reshape(p,N,224)';  %转成M行N列的随机矩阵R CR = 87.5%
% R=reshape(p,N,208)';  %转成M行N列的随机矩阵R CR = 81.25%
% R=reshape(p,N,192)';  %转成M行N列的随机矩阵R CR = 75%
% R=reshape(p,N,128)';  %转成M行N列的随机矩阵R CR = 50%
% R=reshape(p,N,104)';  %转成M行N列的随机矩阵R CR = 40%
R=reshape(p,N,72)';  %转成M行N列的随机矩阵R CR = 20%


%% 4.求解分数阶Chen氏超混沌系统
%求四个初值X0,Y0,Z0,H0


% r=(M/t)*(N/t);      %r为分块个数 CR = 100%
% r=(224/t)*(N/t);      %r为分块个数 CR = 87.5%
% r=(M/t)*(208/t);      %r为分块个数 CR = 81.25%
% r=(M/t)*(192/t);      %r为分块个数 CR = 75%
% r=(M/t)*(128/t);      %r为分块个数 CR = 50%
% r=(M/t)*(104/t);      %r为分块个数 CR = 40%
r=(M/t)*(72/t);      %r为分块个数 CR = 20%


%求出四个初值
X0=sum(sum(bitand(I1,17)))/(17*255*256);
Y0=sum(sum(bitand(I2,34)))/(34*255*256);
Z0=sum(sum(bitand(I3,68)))/(68*255*256);
H0=sum(sum(bitand(I1,136)))/(136*255*256);
%保留四位小数
X0=round(X0*10^4)/10^4;%作为密钥6
Y0=round(Y0*10^4)/10^4;%作为密钥7
Z0=round(Z0*10^4)/10^4;%作为密钥8
H0=round(H0*10^4)/10^4;%作为密钥10

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
% H=H(5000:length(H));
H=H(5000:5000+r);

%% 压缩感知处理
U1=I1;
U2=I2;
U3=I3;
F1=im2double(U1);   %把image转换成双精度类型
F2=im2double(U2);   %把image转换成双精度类型
F3=im2double(U3);   %把image转换成双精度类型
%  小波变换矩阵生成
ww=DWT(N);
%  小波变换让图像稀疏化（注意该步骤会耗费时间，但是会增大稀疏度）
B1=ww*sparse(F1)*ww'; %sparse(X)函数：生成稀疏矩阵，将矩阵A转化为稀疏矩阵形式，即矩阵A中任何0元素被去除，非零元素及其下标组成矩阵S。如果A本身是稀疏的，sparse(S)返回S。
B1=full(B1); %把稀疏矩阵X转换为全矩阵存储形式X，还原X

B2=ww*sparse(F2)*ww'; %sparse(X)函数：生成稀疏矩阵，将矩阵A转化为稀疏矩阵形式，即矩阵A中任何0元素被去除，非零元素及其下标组成矩阵S。如果A本身是稀疏的，sparse(S)返回S。
B2=full(B2); %把稀疏矩阵X转换为全矩阵存储形式X，还原X

B3=ww*sparse(F3)*ww'; %sparse(X)函数：生成稀疏矩阵，将矩阵A转化为稀疏矩阵形式，即矩阵A中任何0元素被去除，非零元素及其下标组成矩阵S。如果A本身是稀疏的，sparse(S)返回S。
B3=full(B3); %把稀疏矩阵X转换为全矩阵存储形式X，还原X


%  随机矩阵生成
% M2=256;      % CR = 100%
% M2=224;      % CR = 87.5%
% M2=208;      % CR = 81.25%
% M2=192;      % CR = 75%
% M2=128;      % CR = 50%
% M2=104;      % CR = 40%
M2=72;       % CR = 20%


K = PartHadamardMtx(M2,N);
save('PartHadamardMatrix.mat','K','ww');

%  测量
Y1 = K*B1;
Y2 = K*B2;
Y3 = K*B3;

save('Measure_MTX.mat','Y1','Y2','Y3');

 %% 数据量化和处理
INTY1=fix(Y1);                                       %对Y1向0取整
INTY2=fix(Y2);                                       %对Y2向0取整
INTY3=fix(Y3);                                       %对Y3向0取整
High_frequency1=Y1-INTY1;                            %高频系数矩阵
High_frequency2=Y2-INTY2;                            %高频系数矩阵
High_frequency3=Y3-INTY3;                            %高频系数矩阵
Index_Matrix1=sign(INTY1);                           %求解正负数矩阵
Index_Matrix2=sign(INTY2);                           %求解正负数矩阵
Index_Matrix3=sign(INTY3);                           %求解正负数矩阵
J1=uint8(abs(INTY1));                                %使得J1中每一个元素均在0~255之间
J2=uint8(abs(INTY2));                                %使得J2中每一个元素均在0~255之间
J3=uint8(abs(INTY3));                                %使得J3中每一个元素均在0~255之间
save('Coefficient_Matrix.mat','High_frequency1','High_frequency2','High_frequency3','Index_Matrix1','Index_Matrix2','Index_Matrix3');
%% 5.DNA编码
%X,Y分别决定I和R的DNA编码方式，有8种，1~8
%Z决定运算方式，有4种，0~3，0表示加，1表示减，2表示异或，3表示同或
%H表示DNA解码方式，有8种，1~8
X=mod(round(X*10^4),8)+1;%使混沌序列X变成范围为1~8的整数
Y=mod(round(Y*10^4),8)+1;%使混沌序列Y变成范围为1~8的整数
Z=mod(round(Z*10^4),4);  %使混沌序列Z变成范围为1~4的整数
H=mod(round(H*10^4),8)+1;%使混沌序列H变成范围为1~8的整数
e=N/t;                   %e表示每一行可以分为多少块

Q2=DNA_bian(fenkuai(t,R,1),Y(1));%调用函数，初始化R矩阵第一块的DNA编码
%R通道
% Q1_R=DNA_bian(fenkuai(t,I1,1),X(1));
Q1_R=DNA_bian(fenkuai(t,J1,1),X(1));%这里对转换量化后的彩色空间Y分量图像进行DNA编码，返回ASCii码
Q_last_R=DNA_yunsuan(Q1_R,Q2,Z(1));
Q_R(1:t,1:t)=DNA_jie(Q_last_R,H(1));
%G通道
% Q1_G=DNA_bian(fenkuai(t,I2,1),X(1));
Q1_G=DNA_bian(fenkuai(t,J2,1),X(1));%这里对转换量化后的彩色空间Cb分量图像进行DNA编码
Q_last_G=DNA_yunsuan(Q1_G,Q2,Z(1));
Q_G(1:t,1:t)=DNA_jie(Q_last_G,H(1));
%B通道
% Q1_B=DNA_bian(fenkuai(t,I3,1),X(1));
Q1_B=DNA_bian(fenkuai(t,J3,1),X(1));%这里对转换量化后的彩色空间Cr分量图像进行DNA编码
Q_last_B=DNA_yunsuan(Q1_B,Q2,Z(1));
Q_B(1:t,1:t)=DNA_jie(Q_last_B,H(1));

for i=2:r
    Q1_R=DNA_bian(fenkuai(t,J1,i),X(i));   %对原始图像R通道每一个分块按X对应的序号进行DNA编码
    Q1_G=DNA_bian(fenkuai(t,J2,i),X(i));   %对原始图像G通道每一个分块按X对应的序号进行DNA编码
    Q1_B=DNA_bian(fenkuai(t,J3,i),X(i));   %对原始图像B通道每一个分块按X对应的序号进行DNA编码
    
    Q2=DNA_bian(fenkuai(t,R,i),Y(i));         %对R的每一个分块按Y对应的序号进行DNA编码
    %R通道->Y通道
    Q3_R=DNA_yunsuan(Q1_R,Q2,Z(i));           %对上面两个编码好的块按Z对应的序号进行DNA运算
    Q4_R=DNA_yunsuan(Q3_R,Q_last_R,Z(i));     %运算结果在和前一块按Z对应的序号再一次进行运算，称为扩散
    Q_last_R=Q4_R;
    %G通道->Cb通道
    Q3_G=DNA_yunsuan(Q1_G,Q2,Z(i));
    Q4_G=DNA_yunsuan(Q3_G,Q_last_G,Z(i));
    Q_last_G=Q4_G;
    %B通道->Cr通道
    Q3_B=DNA_yunsuan(Q1_B,Q2,Z(i));
    Q4_B=DNA_yunsuan(Q3_B,Q_last_B,Z(i));
    Q_last_B=Q4_B;
    
    xx=floor(i/e)+1;
    yy=mod(i,e);
    if yy==0
        xx=xx-1;
        yy=e;
    end
    Q_R((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_R,H(i));    %将每一块合并成完整的图Q，H决定DNA解码过程
    Q_G((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_G,H(i));
    Q_B((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_B,H(i));
end
Q_R=uint8(Q_R);%强制保证Q_R中的值在0~255之间
Q_G=uint8(Q_G);
Q_B=uint8(Q_B);
%% 抗裁剪
xx0=sum(I2(:))/(255*255*256);     %G通道：平均灰度值，作为密钥10
xx0=floor(xx0*10^4)/10^4;     %保留4位小数
xx1=sum(I3(:))/(255*255*256);     %B通道：平均灰度值，作为密钥11
xx1=floor(xx1*10^4)/10^4;     %保留4位小数


% ppx=zeros(1,M+1000);        %预分配内存 %CR = 100%
% ppx=zeros(1,224+1000);        %预分配内存 %CR = 87.5%
% ppx=zeros(1,208+1000);        %预分配内存 %CR = 81.25%
% ppx=zeros(1,192+1000);        %预分配内存 %CR = 75%
% ppx=zeros(1,128+1000);        %预分配内存 %CR = 50%
% ppx=zeros(1,104+1000);        %预分配内存 %CR = 40%
ppx=zeros(1,72+1000);         %预分配内存 %CR = 20%


ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;


% for i=1:M+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
% for i=1:224+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
% for i=1:208+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
% for i=1:192+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
% for i=1:128+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
% for i=1:104+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
for i=1:72+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    

    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %去除前1000点，获得更好的随机性
ppy=ppy(1001:length(ppy));

[v,Ux]=sort(ppx,'descend');
[~,Uy]=sort(ppy,'descend');


% for i=1:M %CR = 100%
% for i=1:224 %CR = 87.5%
% for i=1:208 %CR = 81.25%
% for i=1:192 %CR = 75%
% for i=1:128 %CR = 50%
% for i=1:104 %CR = 40%
for i=1:72  %CR = 20%

    
    temp = Q_R(i,:);
    Q_R(i,:) = Q_R(Ux(i),:);
    Q_R(Ux(i),:) = temp;
    temp = Q_G(i,:);
    Q_G(i,:) = Q_G(Ux(i),:);
    Q_G(Ux(i),:) = temp;
    temp = Q_B(i,:);
    Q_B(i,:) = Q_B(Ux(i),:);
    Q_B(Ux(i),:) = temp;
end

for i=1:N
    temp = Q_R(:,i);
    Q_R(:,i) = Q_R(:,Uy(i));
    Q_R(:,Uy(i)) = temp;
    temp = Q_G(:,i);
    Q_G(:,i) = Q_G(:,Uy(i));
    Q_G(:,Uy(i)) = temp;
    temp = Q_B(:,i);
    Q_B(:,i) = Q_B(:,Uy(i));
    Q_B(:,Uy(i)) = temp;
end

L=figure('color',[1 1 1]);
imhist(Q_R);title('加密后R通道直方图');
axis([0 255 0 500]);
frame = getframe(L);
R_jiamizhifangtu=frame2im(frame);
imshow(R_jiamizhifangtu);

C=figure('color',[1 1 1]);
imhist(Q_G);%title('加密后G通道直方图');
axis([0 255 0 500]);
frame = getframe(C);
G_jiamizhifangtu=frame2im(frame);

B=figure('color',[1 1 1]);
imhist(Q_B);%title('加密后B通道直方图');
axis([0 255 0 500]);
frame = getframe(B);
B_jiamizhifangtu=frame2im(frame);

Q_jiami_RGB(:,:,1)=Q_R;
Q_jiami_RGB(:,:,2)=Q_G;
Q_jiami_RGB(:,:,3)=Q_B;

% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后的lena.png','png'); 

% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.875的lena.png','png'); 
% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.8125的lena.png','png');
% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.75的lena.png','png');
% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.5的lena.png','png');
% imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.4的lena.png','png');
imwrite(Q_jiami_RGB,'../原始、加密、解密图片/lena_CS/加密后CR=0.2的lena.png','png');

figure('color',[1 1 1]);imshow(Q_jiami_RGB);title('加密后图片');
KEY=[t,u,x0,X0,Y0,Z0,H0,M1,N1,xx0,xx1];
save('Encryption_Key.mat','KEY');
%% 输出数据信息
disp('加密成功');
disp('密钥：');    
disp(['密钥1：t=',num2str(t),'     密钥2：u=',num2str(u),'    密钥3：x0=',num2str(x0),'    密钥4：x(0)=',num2str(X0),'   密钥5：y(0)=',num2str(Y0) '   密钥6：z(0)=',num2str(Z0)]);
disp(['密钥7：h(0)=',num2str(H0),'   密钥8：M1=',num2str(M1),'   密钥9：N1=',num2str(N1),'   密钥10：xx0=',num2str(xx0),'   密钥11：xx1=',num2str(xx1)]);