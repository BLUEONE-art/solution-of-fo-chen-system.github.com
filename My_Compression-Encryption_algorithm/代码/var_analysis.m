%% 直方图方差分析
%   @author:董昊
%   @date:2020.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
%I=imread('../原始、加密、解密图片/lena/lena.png','png');         %读取图像信息，lena
%I=imread('../原始、加密、解密图片/辣椒/peppers.png','png');       %读取图像信息，辣椒
%I=imread('../原始、加密、解密图片/狒狒/baboon.png','png');         %读取图像信息，狒狒
%I=imread('../原始、加密、解密图片/飞机/airplane.png','png');         %读取图像信息，飞机
%I=imread('../原始、加密、解密图片/4.2.06小船/4.2.06.tiff','tiff');         %读取图像信息
I=imread('../原始、加密、解密图片/汽车/house.tiff','tiff');         %读取图像信息

figure;imshow(I);title('原始图片');
load('Mask.mat');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B

KEY=zeros(1,15);
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
SUM=M*N;

%% 统计加密前R、G、B三个通道的直方图各个灰度值出现频率
Before_Gray_value_frequency_R = zeros(1,256);
Before_Gray_value_frequency_R(:,1) = length(find(I1 == 0));
for a = 1:255
    Before_Gray_value_frequency_R(:,a+1) = length(find(I1 == a));
end

Before_Gray_value_frequency_G = zeros(1,256);
Before_Gray_value_frequency_G(:,1) = length(find(I2 == 0));
for b = 1:255
    Before_Gray_value_frequency_G(:,b+1) = length(find(I2 == b));
end

Before_Gray_value_frequency_B = zeros(1,256);
Before_Gray_value_frequency_B(:,1) = length(find(I3 == 0));
for c = 1:255
    Before_Gray_value_frequency_B(:,c+1) = length(find(I3 == c));
end

%% 计算加密前直方图方差
Before_var_R = 0;
Before_var_G = 0;
Before_var_B = 0;
for m=1:256
    for n=1:256
        Before_var_R=Before_var_R+(Before_Gray_value_frequency_R(:,m)-Before_Gray_value_frequency_R(:,n))^2;       %R通道MSE
        Before_var_G=Before_var_G+(Before_Gray_value_frequency_G(:,m)-Before_Gray_value_frequency_G(:,n))^2;       %G通道MSE
        Before_var_B=Before_var_B+(Before_Gray_value_frequency_B(:,m)-Before_Gray_value_frequency_B(:,n))^2;       %B通道MSE
    end
end

Before_var_R=Before_var_R/SUM;
Before_var_G=Before_var_G/SUM;
Before_var_B=Before_var_B/SUM;

%% 2.产生Logistic混沌序列
u=3.9999;     %Logistic参数μ，自定为3.9999，作为密钥4
x0=(sum(I1(:))+sum(I2(:)))/(255*SUM*2);     %计算得出Logistic初值x0，作为密钥5，这里是根据原始RGB彩色空间图像进行算的密钥，只跟明文相关
% x0=floor(x0*10^10)/10^10;     %保留10位小数
x0=roundn(x0,-10);
p=zeros(1,SUM+1000);            %预分配内存
p(1)=x0;
for i=1:SUM+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));            %去除前1000点，获得更好的随机性

%% 3.将p序列变换到0~255范围内整数，转换成M*N的二维矩阵R,这里要根据压缩率进行修改！
p=mod(round(p*10^4),256);
R=reshape(p,N,M)';  %转成M行N列的随机矩阵R
%% 4.求解分数阶Chen氏超混沌系统
%求四个初值X0,Y0,Z0,H0
r=(M/t)*(N/t);      %r为分块个数
%求出四个初值
X0=sum(sum(bitand(I1,17)))/(17*SUM);
Y0=sum(sum(bitand(I2,34)))/(34*SUM);
Z0=sum(sum(bitand(I3,68)))/(68*SUM);
H0=sum(sum(bitand(I1,136)))/(136*SUM);
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
%% DCT处理
U1=I1;
U2=I2;
U3=I3;
T=dctmtx(t);        %4阶DCT变换矩阵
F1=im2double(U1);   %把image转换成双精度类型
F2=im2double(U2);   %把image转换成双精度类型
F3=im2double(U3);   %把image转换成双精度类型
fun = @(block_struct) T*block_struct.data*T';           %Matlab建议用blockproc函数代替blkproc函数
J1 = blockproc(F1,[t t],fun);                           %得到Y通道的DCT矩阵J1
J2 = blockproc(F2,[t t],fun);                           %得到Cb通道的DCT矩阵J2
J3 = blockproc(F3,[t t],fun);                           %得到Cr通道的DCT矩阵J3.三个矩阵均为实矩阵

%% 数据量化和处理
max_L=zeros(1,3);
max_firt_num_RGB=zeros(1,3);
BLOCK=ones(t,t);
BLOCK(1,1)=0;
fun = @(block_struct) block_struct.data.*BLOCK; %把每个小块的第一行第一列元素去掉

L1=blockproc(J1,[t t],fun);
L2=blockproc(J2,[t t],fun);
L3=blockproc(J3,[t t],fun);

max_L(1,1)=max(max(abs(L1)));                   %求L1中除去每个分块第一个元素后的最大值
max_L(1,2)=max(max(abs(L2)));
max_L(1,3)=max(max(abs(L3)));

max_num=max(max(max_L));                           %求三个通道中除去每个分块第一个元素后的最大值
max_Multiple=fix(255/max_num);                     %求最大的倍数

max_firt_num_RGB(1,1)=max(max(abs(J1))); 
max_firt_num_RGB(1,2)=max(max(abs(J2))); 
max_firt_num_RGB(1,3)=max(max(abs(J3))); 
max_firt_num=max(max(max_firt_num_RGB));
max_Multiple_first=fix(255/max_firt_num);

rate=4;         %压缩矩阵1存在的行列数，可作为密钥
T1=zeros(t,t);  % 压缩矩阵
T1(1,1)=max_Multiple_first;%可设置成密钥
for k=2:rate
    T1(1,k)=max_Multiple;%可设置成密钥
end
for i=2:rate
    for j=1:rate%+1-i
        T1(i,j)=max_Multiple;
    end
end

fun = @(block_struct) block_struct.data.*T1;         %Matlab建议用blockproc函数代替blkproc函数
J1 = blockproc(J1,[t t],fun);                        %经过数据处理之后的矩阵
J2 = blockproc(J2,[t t],fun);                        %经过数据处理之后的矩阵
J3 = blockproc(J3,[t t],fun);                        %经过数据处理之后的矩阵
INTJ1=fix(J1);                                       %对J1向0取整
INTJ2=fix(J2);                                       %对J2向0取整
INTJ3=fix(J3);                                       %对J3向0取整
High_frequency1=J1-INTJ1;                            %高频系数矩阵
High_frequency2=J2-INTJ2;                            %高频系数矩阵
High_frequency3=J3-INTJ3;                            %高频系数矩阵
Index_Matrix1=sign(INTJ1);                           %求解正负数矩阵
Index_Matrix2=sign(INTJ2);                           %求解正负数矩阵
Index_Matrix3=sign(INTJ3);                           %求解正负数矩阵
J1=uint8(abs(INTJ1));                                %使得J1中每一个元素均在0~255之间
J2=uint8(abs(INTJ2));                                %使得J2中每一个元素均在0~255之间
J3=uint8(abs(INTJ3));                                %使得J3中每一个元素均在0~255之间
save('Coefficient_Matrix.mat','High_frequency1','High_frequency2','High_frequency3','Index_Matrix1','Index_Matrix2','Index_Matrix3');
% 5.DNA编码
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
xx0=sum(I2(:))/(255*SUM);     %G通道：平均灰度值，作为密钥10
xx0=floor(xx0*10^4)/10^4;     %保留4位小数
xx1=sum(I3(:))/(255*SUM);     %B通道：平均灰度值，作为密钥11
xx1=floor(xx1*10^4)/10^4;     %保留4位小数
ppx=zeros(1,M+1000);        %预分配内存
ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;
for i=1:M+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %去除前1000点，获得更好的随机性
ppy=ppy(1001:length(ppy));

[v,Ux]=sort(ppx,'descend');
[~,Uy]=sort(ppy,'descend');

for i=1:M
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
end
%% 统计加密后R、G、B三个通道的直方图各个灰度值出现频率
Gray_value_frequency_R = zeros(1,256);
Gray_value_frequency_R(:,1) = length(find(Q_R == 0));
for a = 1:255
    Gray_value_frequency_R(:,a+1) = length(find(Q_R == a));
end

Gray_value_frequency_G = zeros(1,256);
Gray_value_frequency_G(:,1) = length(find(Q_G == 0));
for a = 1:255
    Gray_value_frequency_G(:,a+1) = length(find(Q_G == a));
end

Gray_value_frequency_B = zeros(1,256);
Gray_value_frequency_B(:,1) = length(find(Q_B == 0));
for a = 1:255
    Gray_value_frequency_B(:,a+1) = length(find(Q_B == a));
end

%% 计算加密前直方图方差
After_var_R = 0;
After_var_G = 0;
After_var_B = 0;
for m=1:256
    for n=1:256
        After_var_R=After_var_R+(Gray_value_frequency_R(:,m)-Gray_value_frequency_R(:,n))^2;       %R通道MSE
        After_var_G=After_var_G+(Gray_value_frequency_G(:,m)-Gray_value_frequency_G(:,n))^2;       %G通道MSE
        After_var_B=After_var_B+(Gray_value_frequency_B(:,m)-Gray_value_frequency_B(:,n))^2;       %B通道MSE
    end
end

After_var_R=After_var_R/SUM;
After_var_G=After_var_G/SUM;
After_var_B=After_var_B/SUM;

disp(['加密前R通道直方图方差：',num2str(Before_var_R)]);
disp(['加密后R通道直方图方差：', num2str(After_var_R)]);

disp(['加密前G通道直方图方差：',num2str(Before_var_G)]);
disp(['加密后G通道直方图方差：', num2str(After_var_G)]);

disp(['加密前B通道直方图方差：',num2str(Before_var_B)]);
disp(['加密后B通道直方图方差：', num2str(After_var_B)]);
