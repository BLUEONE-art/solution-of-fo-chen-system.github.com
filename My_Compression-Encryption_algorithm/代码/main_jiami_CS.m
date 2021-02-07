%% ���ϲ�ɫ����ͼ�����-ѹ��ϵͳ������CS���ܣ�
%   @author:���
%   @date:2020.08.07
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');         %��ȡͼ����Ϣ��lena
figure;imshow(I);title('ԭʼͼƬ');
load('Mask.mat');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B

KEY=zeros(1,11);

Q=figure('color',[1 1 1]);
imhist(I1);%title('ԭʼͼƬRͨ��ֱ��ͼ');
frame = getframe(Q);
R_zhifangtu=frame2im(frame);

W=figure('color',[1 1 1]);
imhist(I2);%title('ԭʼͼƬGͨ��ֱ��ͼ');
frame = getframe(W);
G_zhifangtu=frame2im(frame);

E=figure('color',[1 1 1]);
imhist(I3);%title('ԭʼͼƬBͨ��ֱ��ͼ');
frame = getframe(E);
B_zhifangtu=frame2im(frame);

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


% SUM=M*N; %CR = 100%
% SUM=224*N; %CR = 87.5%
% SUM=M*208; %CR = 81.25%
% SUM=M*192; %CR = 75%
% SUM=M*128; %CR = 50%
% SUM=M*104; %CR = 40%
SUM=M*72; %CR = 20%


%% 2.����Logistic��������
u=3.9999;     %Logistic�����̣��Զ�Ϊ3.9999����Ϊ��Կ4
x0=(sum(I1(:))+sum(I2(:)))/(255*256*256*2);     %����ó�Logistic��ֵx0����Ϊ��Կ5�������Ǹ���ԭʼRGB��ɫ�ռ�ͼ����������Կ��ֻ���������
% x0=floor(x0*10^10)/10^10;     %����10λС��
x0=roundn(x0,-10);


% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 100%
% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 87.5%
% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 81.25%
% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 75%
% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 50%
% p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 40%
p=zeros(1,SUM+1000);            %Ԥ�����ڴ� CR = 20%


p(1)=x0;
for i=1:SUM+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));            %ȥ��ǰ1000�㣬��ø��õ������

%% 3.��p���б任��0~255��Χ��������ת����M*N�Ķ�ά����R,����Ҫ����ѹ���ʽ����޸ģ�
p=mod(round(p*10^4),256); 


% R=reshape(p,N,M)';  %ת��M��N�е��������R CR = 100%
% R=reshape(p,N,224)';  %ת��M��N�е��������R CR = 87.5%
% R=reshape(p,N,208)';  %ת��M��N�е��������R CR = 81.25%
% R=reshape(p,N,192)';  %ת��M��N�е��������R CR = 75%
% R=reshape(p,N,128)';  %ת��M��N�е��������R CR = 50%
% R=reshape(p,N,104)';  %ת��M��N�е��������R CR = 40%
R=reshape(p,N,72)';  %ת��M��N�е��������R CR = 20%


%% 4.��������Chen�ϳ�����ϵͳ
%���ĸ���ֵX0,Y0,Z0,H0


% r=(M/t)*(N/t);      %rΪ�ֿ���� CR = 100%
% r=(224/t)*(N/t);      %rΪ�ֿ���� CR = 87.5%
% r=(M/t)*(208/t);      %rΪ�ֿ���� CR = 81.25%
% r=(M/t)*(192/t);      %rΪ�ֿ���� CR = 75%
% r=(M/t)*(128/t);      %rΪ�ֿ���� CR = 50%
% r=(M/t)*(104/t);      %rΪ�ֿ���� CR = 40%
r=(M/t)*(72/t);      %rΪ�ֿ���� CR = 20%


%����ĸ���ֵ
X0=sum(sum(bitand(I1,17)))/(17*255*256);
Y0=sum(sum(bitand(I2,34)))/(34*255*256);
Z0=sum(sum(bitand(I3,68)))/(68*255*256);
H0=sum(sum(bitand(I1,136)))/(136*255*256);
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

%% ѹ����֪����
U1=I1;
U2=I2;
U3=I3;
F1=im2double(U1);   %��imageת����˫��������
F2=im2double(U2);   %��imageת����˫��������
F3=im2double(U3);   %��imageת����˫��������
%  С���任��������
ww=DWT(N);
%  С���任��ͼ��ϡ�軯��ע��ò����ķ�ʱ�䣬���ǻ�����ϡ��ȣ�
B1=ww*sparse(F1)*ww'; %sparse(X)����������ϡ����󣬽�����Aת��Ϊϡ�������ʽ��������A���κ�0Ԫ�ر�ȥ��������Ԫ�ؼ����±���ɾ���S�����A������ϡ��ģ�sparse(S)����S��
B1=full(B1); %��ϡ�����Xת��Ϊȫ����洢��ʽX����ԭX

B2=ww*sparse(F2)*ww'; %sparse(X)����������ϡ����󣬽�����Aת��Ϊϡ�������ʽ��������A���κ�0Ԫ�ر�ȥ��������Ԫ�ؼ����±���ɾ���S�����A������ϡ��ģ�sparse(S)����S��
B2=full(B2); %��ϡ�����Xת��Ϊȫ����洢��ʽX����ԭX

B3=ww*sparse(F3)*ww'; %sparse(X)����������ϡ����󣬽�����Aת��Ϊϡ�������ʽ��������A���κ�0Ԫ�ر�ȥ��������Ԫ�ؼ����±���ɾ���S�����A������ϡ��ģ�sparse(S)����S��
B3=full(B3); %��ϡ�����Xת��Ϊȫ����洢��ʽX����ԭX


%  �����������
% M2=256;      % CR = 100%
% M2=224;      % CR = 87.5%
% M2=208;      % CR = 81.25%
% M2=192;      % CR = 75%
% M2=128;      % CR = 50%
% M2=104;      % CR = 40%
M2=72;       % CR = 20%


K = PartHadamardMtx(M2,N);
save('PartHadamardMatrix.mat','K','ww');

%  ����
Y1 = K*B1;
Y2 = K*B2;
Y3 = K*B3;

save('Measure_MTX.mat','Y1','Y2','Y3');

 %% ���������ʹ���
INTY1=fix(Y1);                                       %��Y1��0ȡ��
INTY2=fix(Y2);                                       %��Y2��0ȡ��
INTY3=fix(Y3);                                       %��Y3��0ȡ��
High_frequency1=Y1-INTY1;                            %��Ƶϵ������
High_frequency2=Y2-INTY2;                            %��Ƶϵ������
High_frequency3=Y3-INTY3;                            %��Ƶϵ������
Index_Matrix1=sign(INTY1);                           %�������������
Index_Matrix2=sign(INTY2);                           %�������������
Index_Matrix3=sign(INTY3);                           %�������������
J1=uint8(abs(INTY1));                                %ʹ��J1��ÿһ��Ԫ�ؾ���0~255֮��
J2=uint8(abs(INTY2));                                %ʹ��J2��ÿһ��Ԫ�ؾ���0~255֮��
J3=uint8(abs(INTY3));                                %ʹ��J3��ÿһ��Ԫ�ؾ���0~255֮��
save('Coefficient_Matrix.mat','High_frequency1','High_frequency2','High_frequency3','Index_Matrix1','Index_Matrix2','Index_Matrix3');
%% 5.DNA����
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
xx0=sum(I2(:))/(255*255*256);     %Gͨ����ƽ���Ҷ�ֵ����Ϊ��Կ10
xx0=floor(xx0*10^4)/10^4;     %����4λС��
xx1=sum(I3(:))/(255*255*256);     %Bͨ����ƽ���Ҷ�ֵ����Ϊ��Կ11
xx1=floor(xx1*10^4)/10^4;     %����4λС��


% ppx=zeros(1,M+1000);        %Ԥ�����ڴ� %CR = 100%
% ppx=zeros(1,224+1000);        %Ԥ�����ڴ� %CR = 87.5%
% ppx=zeros(1,208+1000);        %Ԥ�����ڴ� %CR = 81.25%
% ppx=zeros(1,192+1000);        %Ԥ�����ڴ� %CR = 75%
% ppx=zeros(1,128+1000);        %Ԥ�����ڴ� %CR = 50%
% ppx=zeros(1,104+1000);        %Ԥ�����ڴ� %CR = 40%
ppx=zeros(1,72+1000);         %Ԥ�����ڴ� %CR = 20%


ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;


% for i=1:M+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
% for i=1:224+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
% for i=1:208+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
% for i=1:192+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
% for i=1:128+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
% for i=1:104+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
for i=1:72+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    

    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %ȥ��ǰ1000�㣬��ø��õ������
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
imhist(Q_R);title('���ܺ�Rͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(L);
R_jiamizhifangtu=frame2im(frame);
imshow(R_jiamizhifangtu);

C=figure('color',[1 1 1]);
imhist(Q_G);%title('���ܺ�Gͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(C);
G_jiamizhifangtu=frame2im(frame);

B=figure('color',[1 1 1]);
imhist(Q_B);%title('���ܺ�Bͨ��ֱ��ͼ');
axis([0 255 0 500]);
frame = getframe(B);
B_jiamizhifangtu=frame2im(frame);

Q_jiami_RGB(:,:,1)=Q_R;
Q_jiami_RGB(:,:,2)=Q_G;
Q_jiami_RGB(:,:,3)=Q_B;

% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ��lena.png','png'); 

% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.875��lena.png','png'); 
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.8125��lena.png','png');
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.75��lena.png','png');
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.5��lena.png','png');
% imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.4��lena.png','png');
imwrite(Q_jiami_RGB,'../ԭʼ�����ܡ�����ͼƬ/lena_CS/���ܺ�CR=0.2��lena.png','png');

figure('color',[1 1 1]);imshow(Q_jiami_RGB);title('���ܺ�ͼƬ');
KEY=[t,u,x0,X0,Y0,Z0,H0,M1,N1,xx0,xx1];
save('Encryption_Key.mat','KEY');
%% ���������Ϣ
disp('���ܳɹ�');
disp('��Կ��');    
disp(['��Կ1��t=',num2str(t),'     ��Կ2��u=',num2str(u),'    ��Կ3��x0=',num2str(x0),'    ��Կ4��x(0)=',num2str(X0),'   ��Կ5��y(0)=',num2str(Y0) '   ��Կ6��z(0)=',num2str(Z0)]);
disp(['��Կ7��h(0)=',num2str(H0),'   ��Կ8��M1=',num2str(M1),'   ��Կ9��N1=',num2str(N1),'   ��Կ10��xx0=',num2str(xx0),'   ��Կ11��xx1=',num2str(xx1)]);