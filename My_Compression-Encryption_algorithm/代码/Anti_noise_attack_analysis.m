%% ��������������
%   @author:���
%   @date:2020.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/���ܺ��lena.png','png');           %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��peppers.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��baboon.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/���ܺ��airplane.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/���ܺ��4.1.05.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/���ܺ��4.2.06.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��house.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/���ܺ��kodim08.png','png');             %��ȡͼ����Ϣ
load('Encryption_Key.mat');
load('Mask.mat');
U1=I(:,:,1);  
U2=I(:,:,2);  
U3=I(:,:,3);  

[M,N]=size(U1);                      %��ͼ������и�ֵ��M,N
t=KEY(1,1);    %�ֿ��С
SUM=M*N;
% u=3.9999;
u=KEY(1,2);
xx0=KEY(1,10);
xx1=KEY(1,11);
ppx=zeros(1,M+1000);        %Ԥ�����ڴ�
ppy=zeros(1,N+1000); 
ppx(1)=xx0;
ppy(1)=xx1;
for i=1:M+999                 %����M+999��ѭ�������õ�M+1000�㣨������ֵ��
    ppx(i+1)=u*ppx(i)*(1-ppx(i));
end
for i=1:N+999                 %����M+999��ѭ�������õ�M+1000�㣨������ֵ��
    ppy(i+1)=u*ppy(i)*(1-ppy(i));
end
ppx=ppx(1001:length(ppx));            %ȥ��ǰ1000�㣬��ø��õ������
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
%% 2.����Logistic��������
% u=3.990000000000001; %��Կ�����Բ���  10^-15
% u=3.99;%��Կ��Logistic������
% x0=0.7067000000000001; %��Կ�����Բ���  10^-16
% x0=0.5475; %��Կ��Logistic��ֵx0
 
x0=KEY(1,3);

p=zeros(1,SUM+1000);
p(1)=x0;
for i=1:SUM+999                        %����SUM+999��ѭ��������SUM+1000������
    p(i+1)=u*p(i)*(1-p(i));
end
p=p(1001:length(p));
%% 3.��p���б任��0~255��Χ��������ת����M*N�Ķ�ά����R
p=mod(round(p*10^4),256);
R=reshape(p,N,M)';  %ת��M��N��

%% 4.�����緽��
% ���ĸ���ֵX0,Y0,Z0,H0
r=(M/t)*(N/t);
X0=KEY(1,4);
Y0=KEY(1,5);
Z0=KEY(1,6);
H0=KEY(1,7);
%%������Chenϵͳ
% A=chen_output(X0,Y0,Z0,H0,r);
% X=A(:,1);
% X=X(3002:length(X));
% Y=A(:,2);
% Y=Y(3002:length(Y));
% Z=A(:,3);
% Z=Z(3002:length(Z));
% H=A(:,4);
% H=H(3002:length(H));

%%�����׳�����Chenϵͳ
h=0.01;%���ò���
NN=25000;%���л���ϵͳ�õ������и���
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

 %X,Y�ֱ����I��R��DNA���뷽ʽ����8�֣�1~8
X=mod(round(X*10^4),8)+1;
Y=mod(round(Y*10^4),8)+1;
Z=mod(round(Z*10^4),4);
Z(Z==0)=4;
Z(Z==1)=0;
Z(Z==4)=1;
H=mod(round(H*10^4),8)+1;
e=N/t;
%% ͼ����������
YY=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');        %��ȡͼ����Ϣ
YY=double(YY);
Y1=YY(:,:,1);        %R
Y2=YY(:,:,2);        %G
Y3=YY(:,:,3);        %B
MSE_R=zeros(1,21);MSE_G=zeros(1,21);MSE_B=zeros(1,21);
j=0;        %�����±�
for i=0:5:5
    P = imnoise(I, 'gaussian', 0, i^2/255^2);  %�����˹������
    imwrite(P,'../�����˹�������ܽ��/lena/��˹��������Ϊ5ʱ�ļ���lena.png','png'); 
% for i=0:0.05:0.05
%     P = imnoise(I,'salt & pepper',i);
    U1=P(:,:,1);     %Rͨ��
    U2=P(:,:,2);     %Gͨ��
    U3=P(:,:,3);     %Bͨ��
    j=j+1;      %�����±�

    %% DNA���루���ܣ�
    for ii=r:-1:2
        Q1_Y=DNA_bian(fenkuai(t,U1,ii),H(ii));
        Q1_Cb=DNA_bian(fenkuai(t,U2,ii),H(ii));
        Q1_Cr=DNA_bian(fenkuai(t,U3,ii),H(ii));

        Q1_last_Y=DNA_bian(fenkuai(t,U1,ii-1),H(ii-1));
        Q1_last_Cb=DNA_bian(fenkuai(t,U2,ii-1),H(ii-1));
        Q1_last_Cr=DNA_bian(fenkuai(t,U3,ii-1),H(ii-1));

        Q2_Y=DNA_yunsuan(Q1_Y,Q1_last_Y,Z(ii));        %��ɢǰ
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
    
   %% �����ݴ���
    rate=4;     %ѹ������1���ڵ�������������Ϊ��Կ
    max_Multiple_first=KEY(1,12);%������Կ֮һ
    max_Multiple=KEY(1,13);%������Կ֮һ
    load('Coefficient_Matrix.mat');
    T1=zeros(t,t);  % ѹ������
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
%     T1=T1.*YS1;%û��ѹ�������
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
    fun = @(block_struct) block_struct.data./T1;         %Matlab������blockproc��������blkproc����
    J1 = blockproc(U1,[t t],fun);                        %�������ݴ���֮��ľ���
    J2 = blockproc(U2,[t t],fun);                        %�������ݴ���֮��ľ���
    J3 = blockproc(U3,[t t],fun);                        %�������ݴ���֮��ľ���
    %% ��DCT����
    T=dctmtx(t);        %8��DCT�任����
    fun = @(block_struct) T'*block_struct.data*T;       %Matlab������blockproc��������blkproc����
    K1 = blockproc(J1,[t t],fun);              %DCT��任���ع�ͼ��
    K2 = blockproc(J2,[t t],fun);              %DCT��任���ع�ͼ��
    K3 = blockproc(J3,[t t],fun);              %DCT��任���ع�ͼ��

%     I_jiemi(:,:,1)=K1;
%     I_jiemi(:,:,2)=K2;
%     I_jiemi(:,:,3)=K3;

    K1=im2uint8(K1);
    K2=im2uint8(K2);
    K3=im2uint8(K3);
    
    K1=double(K1);
    K2=double(K2);
    K3=double(K3);
    
    %% 6��ȥ������ʱ������
    M1=KEY(1,8);   %����ʱ����Ĳ�����M1=mod(M,t);��Ϊ��Կ
    N1=KEY(1,9);   %����ʱ����Ĳ�����N1=mod(N,t);��Ϊ��Կ
% %     if M1~=0
% %         I_jiemi=I_jiemi(1:M-t+M1,:,:);
% %     end
% %     if N1~=0
% %         I_jiemi=I_jiemi(:,1:N-t+N1,:);
% %     end
    %ȥ������ʱ������
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
    [MM,NN]=size(K1);     %���»�ý��ܺ��ͼƬ��С
    for m=1:MM
        for n=1:NN
            MSE_R(j)=MSE_R(j)+(Y1(m,n)-K1(m,n))^2;       %Rͨ��MSE
            MSE_G(j)=MSE_G(j)+(Y2(m,n)-K2(m,n))^2;       %Gͨ��MSE
            MSE_B(j)=MSE_B(j)+(Y3(m,n)-K3(m,n))^2;       %Bͨ��MSE
        end
    end

    I_jiemi(:,:,1)=uint8(K1);
    I_jiemi(:,:,2)=uint8(K2);
    I_jiemi(:,:,3)=uint8(K3);
%     I_jiemi(:,:,1)=K1;
%     I_jiemi(:,:,2)=K2;
%     I_jiemi(:,:,3)=K3;
% 
    figure;imshow(I_jiemi);title(['��˹����ǿ��Ϊ',num2str(i),'ʱ�Ľ���ͼ��']);
    imwrite(I_jiemi,'../�����˹�������ܽ��/lena/��˹��������Ϊ5ʱ��lena����.png','png'); 
%     figure;imshow(I_jiemi);title(['���������ܶ�Ϊ',num2str(i),'ʱ�Ľ���ͼ��']);
%     imwrite(I_jiemi,'../���뽷���������ܽ��/lena/��˹��������Ϊ5ʱ��lena����.png','png'); 
end
%��������-MSE
MSE_R=MSE_R./SUM;
MSE_G=MSE_G./SUM;
MSE_B=MSE_B./SUM;
%��ֵ�����-PSNR
PSNR_R=10*log10((255^2)./MSE_R);
PSNR_G=10*log10((255^2)./MSE_G);
PSNR_B=10*log10((255^2)./MSE_B);
% ��ͼ����������-MSE����ֵ�����-PSNR
% set(gca,'xtick', X);
% X=0:0.05:1;
% X=0:5:40;
% figure;plot(X,MSE_R);xlabel('���������ܶ�');ylabel('�������MSE');title('Rͨ�������������ܶ�-�������MSE����ͼ');
% figure;plot(X,MSE_G);xlabel('���������ܶ�');ylabel('�������MSE');title('Gͨ�������������ܶ�-�������MSE����ͼ');
% figure;plot(X,MSE_B);xlabel('���������ܶ�');ylabel('�������MSE');title('Bͨ�������������ܶ�-�������MSE����ͼ');
% figure;plot(X,PSNR_R);xlabel('���������ܶ�');ylabel('��ֵ�����PSNR��dB��');title('Rͨ�������������ܶ�-��ֵ�����PSNR����ͼ');
% figure;plot(X,PSNR_G);xlabel('���������ܶ�');ylabel('��ֵ�����PSNR��dB��');title('Gͨ�������������ܶ�-��ֵ�����PSNR����ͼ');
% figure;plot(X,PSNR_B);xlabel('���������ܶ�');ylabel('��ֵ�����PSNR��dB��');title('Bͨ�������������ܶ�-��ֵ�����PSNR����ͼ');