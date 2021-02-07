%% ֱ��ͼ�������
%   @author:���
%   @date:2020.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
%I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');         %��ȡͼ����Ϣ��lena
%I=imread('../ԭʼ�����ܡ�����ͼƬ/����/peppers.png','png');       %��ȡͼ����Ϣ������
%I=imread('../ԭʼ�����ܡ�����ͼƬ/����/baboon.png','png');         %��ȡͼ����Ϣ������
%I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/airplane.png','png');         %��ȡͼ����Ϣ���ɻ�
%I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/4.2.06.tiff','tiff');         %��ȡͼ����Ϣ
I=imread('../ԭʼ�����ܡ�����ͼƬ/����/house.tiff','tiff');         %��ȡͼ����Ϣ

figure;imshow(I);title('ԭʼͼƬ');
load('Mask.mat');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B

KEY=zeros(1,15);
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

%% ͳ�Ƽ���ǰR��G��B����ͨ����ֱ��ͼ�����Ҷ�ֵ����Ƶ��
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

%% �������ǰֱ��ͼ����
Before_var_R = 0;
Before_var_G = 0;
Before_var_B = 0;
for m=1:256
    for n=1:256
        Before_var_R=Before_var_R+(Before_Gray_value_frequency_R(:,m)-Before_Gray_value_frequency_R(:,n))^2;       %Rͨ��MSE
        Before_var_G=Before_var_G+(Before_Gray_value_frequency_G(:,m)-Before_Gray_value_frequency_G(:,n))^2;       %Gͨ��MSE
        Before_var_B=Before_var_B+(Before_Gray_value_frequency_B(:,m)-Before_Gray_value_frequency_B(:,n))^2;       %Bͨ��MSE
    end
end

Before_var_R=Before_var_R/SUM;
Before_var_G=Before_var_G/SUM;
Before_var_B=Before_var_B/SUM;

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
end
%% ͳ�Ƽ��ܺ�R��G��B����ͨ����ֱ��ͼ�����Ҷ�ֵ����Ƶ��
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

%% �������ǰֱ��ͼ����
After_var_R = 0;
After_var_G = 0;
After_var_B = 0;
for m=1:256
    for n=1:256
        After_var_R=After_var_R+(Gray_value_frequency_R(:,m)-Gray_value_frequency_R(:,n))^2;       %Rͨ��MSE
        After_var_G=After_var_G+(Gray_value_frequency_G(:,m)-Gray_value_frequency_G(:,n))^2;       %Gͨ��MSE
        After_var_B=After_var_B+(Gray_value_frequency_B(:,m)-Gray_value_frequency_B(:,n))^2;       %Bͨ��MSE
    end
end

After_var_R=After_var_R/SUM;
After_var_G=After_var_G/SUM;
After_var_B=After_var_B/SUM;

disp(['����ǰRͨ��ֱ��ͼ���',num2str(Before_var_R)]);
disp(['���ܺ�Rͨ��ֱ��ͼ���', num2str(After_var_R)]);

disp(['����ǰGͨ��ֱ��ͼ���',num2str(Before_var_G)]);
disp(['���ܺ�Gͨ��ֱ��ͼ���', num2str(After_var_G)]);

disp(['����ǰBͨ��ֱ��ͼ���',num2str(Before_var_B)]);
disp(['���ܺ�Bͨ��ֱ��ͼ���', num2str(After_var_B)]);
