%% ���ϲ�ɫ����ͼ�����-ѹ��ϵͳ�����ܣ�
%   @author:���
%   @date:2020.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');         %��ȡͼ����Ϣ��lena
% I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/ԭʼlena/lena.png','png');         %��ȡͼ����Ϣ��lena


% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/peppers.png','png');       %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/baboon.png','png');         %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/airplane.png','png');         %��ȡͼ����Ϣ���ɻ�
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/4.1.05.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/4.2.06.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/house.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/kodim08.png','png');         %��ȡͼ����Ϣ
figure;imshow(I);title('ԭʼͼƬ');
load('Mask.mat');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B

KEY=zeros(1,15);

Q=figure('color',[1 1 1]);
imhist(I1);%title('ԭʼͼƬRͨ��ֱ��ͼ');
frame = getframe(Q);
R_zhifangtu=frame2im(frame);
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/lena/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰRͨ��ֱ��ͼ.png','png');
% imwrite(R_zhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/����ǰRͨ��ֱ��ͼ.png','png');

W=figure('color',[1 1 1]);
imhist(I2);%title('ԭʼͼƬGͨ��ֱ��ͼ');
frame = getframe(W);
G_zhifangtu=frame2im(frame);
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/lena/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰGͨ��ֱ��ͼ.png','png');
% imwrite(G_zhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/����ǰGͨ��ֱ��ͼ.png','png');

E=figure('color',[1 1 1]);
imhist(I3);%title('ԭʼͼƬBͨ��ֱ��ͼ');
frame = getframe(E);
B_zhifangtu=frame2im(frame);
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/lena/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/����/����ǰBͨ��ֱ��ͼ.png','png');
% imwrite(B_zhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/����ǰBͨ��ֱ��ͼ.png','png');

% %%�ı�һ������ֵ������һ������ͼ���Լ�⿹��ֹ�������
% I1(1,1)=238;
% I2(1,1)=238;
% I3(1,1)=238;

[M,N]=size(I1);                      %��ͼ������и�ֵ��M,N
t=4;    %�ֿ��С,��Կ1
%% 1.����
%��ͼ��������������ɿ��Ա�t����������tΪ�ֿ�Ĵ�С��
M1=mod(M,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0����Կ2
N1=mod(N,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0����Կ3
                %����M=5����M1Ϊ1���������油�����
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
[M,N]=size(I1);  %����������������
SUM=M*N;
%% 2.����Logistic��������
u=3.9999;     %Logistic�����̣��Զ�Ϊ3.9999����Ϊ��Կ4
x0=(sum(I1(:))+sum(I2(:)))/(255*SUM*2);     %����ó�Logistic��ֵx0����Ϊ��Կ5�������Ǹ���ԭʼRGB��ɫ�ռ�ͼ����������Կ��ֻ���������
% x0=floor(x0*10^10)/10^10;     %����10λС��
x0=roundn(x0,-10);
p=zeros(1,SUM+1000);            %Ԥ�����ڴ�
p(1)=x0;
for i=1:SUM+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));            %ȥ��ǰ1000�㣬��ø��õ������

%% 3.��p���б任��0~255��Χ��������ת����M*N�Ķ�ά����R,����Ҫ����ѹ���ʽ����޸ģ�
p=mod(round(p*10^4),256);
R=reshape(p,N,M)';  %ת��M��N�е��������R
% fun = @(block_struct) block_struct.data.*YS1;        %ѹ����0.875
% fun = @(block_struct) block_struct.data.*YS2;        %ѹ����0.8125
% fun = @(block_struct) block_struct.data.*YS3;        %ѹ����0.75
% fun = @(block_struct) block_struct.data.*YS4;        %ѹ����0.5
% fun = @(block_struct) block_struct.data.*YS5;        %ѹ����0.4
% fun = @(block_struct) block_struct.data.*YS6;        %ѹ����0.2
% R = blockproc(R,[t t],fun);                        %����ѹ��֮��ľ���
%% 4.��������Chen�ϳ�����ϵͳ
%���ĸ���ֵX0,Y0,Z0,H0
r=(M/t)*(N/t);      %rΪ�ֿ����
%����ĸ���ֵ
X0=sum(sum(bitand(I1,17)))/(17*SUM);
Y0=sum(sum(bitand(I2,34)))/(34*SUM);
Z0=sum(sum(bitand(I3,68)))/(68*SUM);
H0=sum(sum(bitand(I1,136)))/(136*SUM);
%������λС��
X0=round(X0*10^4)/10^4;%��Ϊ��Կ6
Y0=round(Y0*10^4)/10^4;%��Ϊ��Կ7
Z0=round(Z0*10^4)/10^4;%��Ϊ��Կ8
H0=round(H0*10^4)/10^4;%��Ϊ��Կ10
%���ݳ�ֵ�����Chen�ϳ�����ϵͳ���õ��ĸ���������
% h=0.01;%���ò���
% NN=25000;%���л���ϵͳ�õ������и���
% q=0.95;%�����׽���
% z0=[X0 Y0 Z0 H0];
% [~,y]=FrataSim(h,NN,z0,q);
% X=y(1,:);
% X=X(5000:5000+r);        %ȥ��ǰ3001���ø��õ�����ԣ�������ϵͳ���Ӻ����������3000�㣩
% Y=y(2,:);
% Y=Y(5000:5000+r);
% Z=y(3,:);
% Z=Z(5000:5000+r);
% H=y(4,:);
% % H=H(5000:length(H));
% H=H(5000:5000+r);

%��Գ���512*512�ĸ����ͼ����Ҫ����Ļ�������
h=0.01;%���ò���
NN=35000;%���л���ϵͳ�õ������и���
q=0.95;%�����׽���
z0=[X0 Y0 Z0 H0];
[~,y]=FrataSim(h,NN,z0,q);
X=y(1,:);
X=X(5000:5000+r);        %ȥ��ǰ3001���ø��õ�����ԣ�������ϵͳ���Ӻ����������3000�㣩
Y=y(2,:);
Y=Y(5000:5000+r);
Z=y(3,:);
Z=Z(5000:5000+r);
H=y(4,:);
% H=H(5000:length(H));
H=H(5000:5000+r);

% A=chen_output(X0,Y0,Z0,H0,r);   %���������׳�����Chenϵͳ
% X=A(:,1);
% X=X(3002:length(X));        %ȥ��ǰ3001���ø��õ�����ԣ�������ϵͳ���Ӻ����������3000�㣩
% Y=A(:,2);
% Y=Y(3002:length(Y));
% Z=A(:,3);
% Z=Z(3002:length(Z));
% H=A(:,4);
% % H=H(5000:length(H));
% H=H(3002:length(H));
%% DCT����
U1=I1;
U2=I2;
U3=I3;
T=dctmtx(t);        %4��DCT�任����
F1=im2double(U1);   %��imageת����˫��������
F2=im2double(U2);   %��imageת����˫��������
F3=im2double(U3);   %��imageת����˫��������
fun = @(block_struct) T*block_struct.data*T';           %Matlab������blockproc��������blkproc����
J1 = blockproc(F1,[t t],fun);                           %�õ�Yͨ����DCT����J1
J2 = blockproc(F2,[t t],fun);                           %�õ�Cbͨ����DCT����J2
J3 = blockproc(F3,[t t],fun);                           %�õ�Crͨ����DCT����J3.���������Ϊʵ����

%% ��������ȫ�������൱�ڲ�������������Ӧ��CR=1
%YS1:CR=0.875
%YS2:CR=0.8125
%YS3:CR=0.75
%YS4:CR=0.5
%YS5:CR=0.4
%YS6:CR=0.2

% a = sum(YS1(:));     %ÿ���ӿ鱣�������ݸ���
% a = sum(YS2(:));     %ÿ���ӿ鱣�������ݸ���
% a = sum(YS3(:));     %ÿ���ӿ鱣�������ݸ���
% a = sum(YS4(:));     %ÿ���ӿ鱣�������ݸ���
% a = sum(YS5(:));     %ÿ���ӿ鱣�������ݸ���
% a = sum(YS6(:));     %ÿ���ӿ鱣�������ݸ���

% fun = @(block_struct) block_struct.data.*YS1;      
% fun = @(block_struct) block_struct.data.*YS2;     
% fun = @(block_struct) block_struct.data.*YS3;    
% fun = @(block_struct) block_struct.data.*YS4;  
% fun = @(block_struct) block_struct.data.*YS5; 
% fun = @(block_struct) block_struct.data.*YS6; 

% J1 = blockproc(J1,[t t],fun);                        %����ѹ��֮��ľ���
% J2 = blockproc(J2,[t t],fun);                        %����ѹ��֮��ľ���
% J3 = blockproc(J3,[t t],fun);                        %����ѹ��֮��ľ���

%% ���������ʹ���
max_L=zeros(1,3);
max_firt_num_RGB=zeros(1,3);
BLOCK=ones(t,t);
BLOCK(1,1)=0;
fun = @(block_struct) block_struct.data.*BLOCK; %��ÿ��С��ĵ�һ�е�һ��Ԫ��ȥ��

L1=blockproc(J1,[t t],fun);
L2=blockproc(J2,[t t],fun);
L3=blockproc(J3,[t t],fun);

max_L(1,1)=max(max(abs(L1)));                   %��L1�г�ȥÿ���ֿ��һ��Ԫ�غ�����ֵ
max_L(1,2)=max(max(abs(L2)));
max_L(1,3)=max(max(abs(L3)));

max_num=max(max(max_L));                           %������ͨ���г�ȥÿ���ֿ��һ��Ԫ�غ�����ֵ
max_Multiple=fix(255/max_num);                     %�����ı���

max_firt_num_RGB(1,1)=max(max(abs(J1))); 
max_firt_num_RGB(1,2)=max(max(abs(J2))); 
max_firt_num_RGB(1,3)=max(max(abs(J3))); 
max_firt_num=max(max(max_firt_num_RGB));
max_Multiple_first=fix(255/max_firt_num);

rate=4;         %ѹ������1���ڵ�������������Ϊ��Կ
T1=zeros(t,t);  % ѹ������
T1(1,1)=max_Multiple_first;%�����ó���Կ
for k=2:rate
    T1(1,k)=max_Multiple;%�����ó���Կ
end
for i=2:rate
    for j=1:rate%+1-i
        T1(i,j)=max_Multiple;
    end
end

% T1=T1.*YS1;
% T1=T1.*YS2;
% T1=T1.*YS3;
% T1=T1.*YS4;
% T1=T1.*YS5;
% T1=T1.*YS6;

fun = @(block_struct) block_struct.data.*T1;         %Matlab������blockproc��������blkproc����
J1 = blockproc(J1,[t t],fun);                        %�������ݴ���֮��ľ���
J2 = blockproc(J2,[t t],fun);                        %�������ݴ���֮��ľ���
J3 = blockproc(J3,[t t],fun);                        %�������ݴ���֮��ľ���
INTJ1=fix(J1);                                       %��J1��0ȡ��
INTJ2=fix(J2);                                       %��J2��0ȡ��
INTJ3=fix(J3);                                       %��J3��0ȡ��
High_frequency1=J1-INTJ1;                            %��Ƶϵ������
High_frequency2=J2-INTJ2;                            %��Ƶϵ������
High_frequency3=J3-INTJ3;                            %��Ƶϵ������
Index_Matrix1=sign(INTJ1);                           %�������������
Index_Matrix2=sign(INTJ2);                           %�������������
Index_Matrix3=sign(INTJ3);                           %�������������
J1=uint8(abs(INTJ1));                                %ʹ��J1��ÿһ��Ԫ�ؾ���0~255֮��
J2=uint8(abs(INTJ2));                                %ʹ��J2��ÿһ��Ԫ�ؾ���0~255֮��
J3=uint8(abs(INTJ3));                                %ʹ��J3��ÿһ��Ԫ�ؾ���0~255֮��
save('Coefficient_Matrix.mat','High_frequency1','High_frequency2','High_frequency3','Index_Matrix1','Index_Matrix2','Index_Matrix3');
% 5.DNA����
%X,Y�ֱ����I��R��DNA���뷽ʽ����8�֣�1~8
%Z�������㷽ʽ����4�֣�0~3��0��ʾ�ӣ�1��ʾ����2��ʾ���3��ʾͬ��
%H��ʾDNA���뷽ʽ����8�֣�1~8
X=mod(round(X*10^4),8)+1;%ʹ��������X��ɷ�ΧΪ1~8������
Y=mod(round(Y*10^4),8)+1;%ʹ��������Y��ɷ�ΧΪ1~8������
Z=mod(round(Z*10^4),4);  %ʹ��������Z��ɷ�ΧΪ1~4������
H=mod(round(H*10^4),8)+1;%ʹ��������H��ɷ�ΧΪ1~8������
e=N/t;                   %e��ʾÿһ�п��Է�Ϊ���ٿ�

Q2=DNA_bian(fenkuai(t,R,1),Y(1));%���ú�������ʼ��R�����һ���DNA����
%Rͨ��
% Q1_R=DNA_bian(fenkuai(t,I1,1),X(1));
Q1_R=DNA_bian(fenkuai(t,J1,1),X(1));%�����ת��������Ĳ�ɫ�ռ�Y����ͼ�����DNA���룬����ASCii��
Q_last_R=DNA_yunsuan(Q1_R,Q2,Z(1));
Q_R(1:t,1:t)=DNA_jie(Q_last_R,H(1));
%Gͨ��
% Q1_G=DNA_bian(fenkuai(t,I2,1),X(1));
Q1_G=DNA_bian(fenkuai(t,J2,1),X(1));%�����ת��������Ĳ�ɫ�ռ�Cb����ͼ�����DNA����
Q_last_G=DNA_yunsuan(Q1_G,Q2,Z(1));
Q_G(1:t,1:t)=DNA_jie(Q_last_G,H(1));
%Bͨ��
% Q1_B=DNA_bian(fenkuai(t,I3,1),X(1));
Q1_B=DNA_bian(fenkuai(t,J3,1),X(1));%�����ת��������Ĳ�ɫ�ռ�Cr����ͼ�����DNA����
Q_last_B=DNA_yunsuan(Q1_B,Q2,Z(1));
Q_B(1:t,1:t)=DNA_jie(Q_last_B,H(1));

for i=2:r
    Q1_R=DNA_bian(fenkuai(t,J1,i),X(i));   %��ԭʼͼ��Rͨ��ÿһ���ֿ鰴X��Ӧ����Ž���DNA����
    Q1_G=DNA_bian(fenkuai(t,J2,i),X(i));   %��ԭʼͼ��Gͨ��ÿһ���ֿ鰴X��Ӧ����Ž���DNA����
    Q1_B=DNA_bian(fenkuai(t,J3,i),X(i));   %��ԭʼͼ��Bͨ��ÿһ���ֿ鰴X��Ӧ����Ž���DNA����
    
    Q2=DNA_bian(fenkuai(t,R,i),Y(i));         %��R��ÿһ���ֿ鰴Y��Ӧ����Ž���DNA����
    %Rͨ��->Yͨ��
    Q3_R=DNA_yunsuan(Q1_R,Q2,Z(i));           %��������������õĿ鰴Z��Ӧ����Ž���DNA����
    Q4_R=DNA_yunsuan(Q3_R,Q_last_R,Z(i));     %�������ں�ǰһ�鰴Z��Ӧ�������һ�ν������㣬��Ϊ��ɢ
    Q_last_R=Q4_R;
    %Gͨ��->Cbͨ��
    Q3_G=DNA_yunsuan(Q1_G,Q2,Z(i));
    Q4_G=DNA_yunsuan(Q3_G,Q_last_G,Z(i));
    Q_last_G=Q4_G;
    %Bͨ��->Crͨ��
    Q3_B=DNA_yunsuan(Q1_B,Q2,Z(i));
    Q4_B=DNA_yunsuan(Q3_B,Q_last_B,Z(i));
    Q_last_B=Q4_B;
    
    xx=floor(i/e)+1;
    yy=mod(i,e);
    if yy==0
        xx=xx-1;
        yy=e;
    end
    Q_R((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_R,H(i));    %��ÿһ��ϲ���������ͼQ��H����DNA�������
    Q_G((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_G,H(i));
    Q_B((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q4_B,H(i));
end
Q_R=uint8(Q_R);%ǿ�Ʊ�֤Q_R�е�ֵ��0~255֮��
Q_G=uint8(Q_G);
Q_B=uint8(Q_B);
%% ���ü�
xx0=sum(I2(:))/(255*SUM);     %Gͨ����ƽ���Ҷ�ֵ����Ϊ��Կ10
xx0=floor(xx0*10^4)/10^4;     %����4λС��
xx1=sum(I3(:))/(255*SUM);     %Bͨ����ƽ���Ҷ�ֵ����Ϊ��Կ11
xx1=floor(xx1*10^4)/10^4;     %����4λС��
ppx=zeros(1,M+1000);        %Ԥ�����ڴ�
ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;
for i=1:M+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %ȥ��ǰ1000�㣬��ø��õ������
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
    Q_B(:,Uy(i)) = temp;
end

L=figure('color',[1 1 1]);
imhist(Q_R);%title('���ܺ�Rͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(L);
R_jiamizhifangtu=frame2im(frame);
imshow(R_jiamizhifangtu);
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/lena/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Rͨ��ֱ��ͼ.png','png');
% imwrite(R_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/���ܺ�Rͨ��ֱ��ͼ.png','png');

C=figure('color',[1 1 1]);
imhist(Q_G);%title('���ܺ�Gͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(C);
G_jiamizhifangtu=frame2im(frame);
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/lena/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Gͨ��ֱ��ͼ.png','png');
% imwrite(G_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/���ܺ�Gͨ��ֱ��ͼ.png','png');

B=figure('color',[1 1 1]);
imhist(Q_B);%title('���ܺ�Bͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(B);
B_jiamizhifangtu=frame2im(frame);
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/lena/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�ɻ�/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.1.05����/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/4.2.06С��/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/����/���ܺ�Bͨ��ֱ��ͼ.png','png');
% imwrite(B_jiamizhifangtu,'../����ǰ��ֱ��ͼ/�Ǳ�/���ܺ�Bͨ��ֱ��ͼ.png','png');

Q_jiami_RGB(:,:,1)=Q_R;
Q_jiami_RGB(:,:,2)=Q_G;
Q_jiami_RGB(:,:,3)=Q_B;

% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena/���ܺ��lena.png','png'); 

% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena/ԭʼlena/���ܺ��lena.png','png'); 

% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena/�ı�һ�����ؼ��ܺ��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/�ı�һ�����ؼ��ܺ��peppers.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/�ı�һ�����ؼ��ܺ��baboon.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/�ı�һ�����ؼ��ܺ��airplane.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/�ı�һ�����ؼ��ܺ��house.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/�ı�һ�����ؼ��ܺ��4.2.06.png','png'); 


% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��peppers.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��baboon.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/���ܺ��airplane.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/4.1.05����/���ܺ��4.1.05.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/���ܺ��4.2.06.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��house.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/���ܺ��kodim08.png','png'); 

% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.875��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.8125��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.75��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.5��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.4��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../PSNR����/lena/���ܺ�CR=0.2��lena.png','png'); 

figure('color',[1 1 1]);imshow(Q_jiami_RGB);title('���ܺ�ͼƬ');
KEY=[t,u,x0,X0,Y0,Z0,H0,M1,N1,xx0,xx1,max_Multiple_first,max_Multiple];
save('Encryption_Key.mat','KEY');
%% ���������Ϣ
disp('���ܳɹ�');
disp('��Կ��');    
disp(['��Կ1��t=',num2str(t),'     ��Կ2��u=',num2str(u),'    ��Կ3��x0=',num2str(x0),'    ��Կ4��x(0)=',num2str(X0),'   ��Կ5��y(0)=',num2str(Y0) '   ��Կ6��z(0)=',num2str(Z0)]);
disp(['��Կ7��h(0)=',num2str(H0),'   ��Կ8��M1=',num2str(M1),'   ��Կ9��N1=',num2str(N1),'   ��Կ10��xx0=',num2str(xx0),'   ��Կ11��xx1=',num2str(xx1) '   ��Կ12��max_Multiple_first=',num2str(max_Multiple_first) '   ��Կ13��max_Multiple=',num2str(max_Multiple)]);