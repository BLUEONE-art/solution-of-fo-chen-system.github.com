%%画分数阶超混沌Chen系统的吸引子图
%%时间：5.11.2020
%%作者：董昊

clear;clc;
I=imread('../原始、加密、解密图片/lena/lena.png','png');         %读取图像信息
% I=imread('../原始、加密、解密图片/辣椒/peppers.png','png');       %读取图像信息，辣椒
% I=imread('../原始、加密、解密图片/狒狒/baboon.png','png');         %读取图像信息，狒狒
% I=imread('../原始、加密、解密图片/飞机/airplane.png','png');         %读取图像信息，飞机
% I=imread('../原始、加密、解密图片/4.1.05房子/4.1.05.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/4.2.06小船/4.2.06.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/汽车/house.tiff','tiff');         %读取图像信息
% I=imread('../原始、加密、解密图片/城堡/kodim08.png','png');         %读取图像信息
figure;imshow(I);title('原始图片');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
[M,N]=size(I1);                      %将图像的行列赋值给M,N
t=4;    %分块大小
%% 1.补零
%将图像的行列数都补成可以被t整除的数，t为分块的大小。
M1=mod(M,t);    %可作为固定密钥，以便解码时可以去除补上的0
N1=mod(N,t);    %可作为固定密钥，以便解码时可以去除补上的0
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
%% 4.求解分数阶Chen氏超混沌系统
%求四个初值X0,Y0,Z0,H0
r=(M/t)*(N/t);      %r为分块个数
%求出四个初值
X0=sum(sum(bitand(I1,17)))/(17*SUM);
Y0=sum(sum(bitand(I2,34)))/(34*SUM);
Z0=sum(sum(bitand(I3,68)))/(68*SUM);
H0=sum(sum(bitand(I1,136)))/(136*SUM);
%保留四位小数
X0=round(X0*10^4)/10^4;%作为密钥3
Y0=round(Y0*10^4)/10^4;%作为密钥4
Z0=round(Z0*10^4)/10^4;%作为密钥5
H0=round(H0*10^4)/10^4;%作为密钥6
%根据初值，求解Chen氏超混沌系统，得到四个混沌序列  
%超混沌Chen系统
% A=chen_output(X0,Y0,Z0,H0,r);
% X=A(:,1);
% X=X(3002:length(X));
% Y=A(:,2);
% Y=Y(3002:length(Y));
% Z=A(:,3);
% Z=Z(3002:length(Z));
% H=A(:,4);
% H=H(3002:length(H));

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

%%混沌吸引子图
Q=figure(1);
plot(X,Y);
xlabel('\itx')
ylabel('\ity')
frame = getframe(Q);
x_y=frame2im(frame);
imwrite(x_y,'../混沌图/lena/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/辣椒/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/狒狒/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/飞机/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/4.1.05房子/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/4.2.06小船/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/汽车/分数阶超混沌Chen系统吸引子x-y.png','png');
% imwrite(x_y,'../混沌图/城堡/分数阶超混沌Chen系统吸引子x-y.png','png');

W=figure(2);
plot(X,Z)
xlabel('\itx')
ylabel('\itz')
frame = getframe(W);
x_z=frame2im(frame);
imwrite(x_z,'../混沌图/lena/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/辣椒/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/狒狒/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/飞机/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/4.1.05房子/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/4.2.06小船/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/汽车/分数阶超混沌Chen系统吸引子x-z.png','png');
% imwrite(x_z,'../混沌图/城堡/分数阶超混沌Chen系统吸引子x-z.png','png');

E=figure(3);
plot(X,H)
xlabel('\itx')
ylabel('\itw')
frame = getframe(E);
x_w=frame2im(frame);
imwrite(x_w,'../混沌图/lena/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/辣椒/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/狒狒/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/飞机/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/4.1.05房子/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/4.2.06小船/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/汽车/分数阶超混沌Chen系统吸引子x-w.png','png');
% imwrite(x_w,'../混沌图/城堡/分数阶超混沌Chen系统吸引子x-w.png','png');

R=figure(4);
plot(Z,H)
xlabel('\itz')
ylabel('\itw')
frame = getframe(R);
z_w=frame2im(frame);
imwrite(z_w,'../混沌图/lena/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/辣椒/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/狒狒/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/飞机/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/4.1.05房子/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/4.2.06小船/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/汽车/分数阶超混沌Chen系统吸引子z-w.png','png');
% imwrite(z_w,'../混沌图/城堡/分数阶超混沌Chen系统吸引子z-w.png','png');