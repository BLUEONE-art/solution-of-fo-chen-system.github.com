%% ���ϲ�ɫ����ͼ�����-ѹ��ϵͳ�����ܣ�
%   @author:���
%   @date:2018.04.17
%-------------------------------------------------------------------------------------------------------%
clear;clc;
I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/���ܺ��lena.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/ԭʼlena/���ܺ��lena.png','png');             %��ȡͼ����Ϣ

% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��peppers.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��baboon.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/���ܺ��airplane.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/���ܺ��4.1.05.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/���ܺ��4.2.06.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��house.png','png');             %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/���ܺ��kodim08.png','png');             %��ȡͼ����Ϣ

% I=imread('../PSNR����/lena/���ܺ�CR=1��lena.png','png'); 
% I=imread('../PSNR����/lena/���ܺ�CR=0.875��lena.png','png'); 
% I=imread('../PSNR����/lena/���ܺ�CR=0.75��lena.png','png'); 
% I=imread('../PSNR����/lena/���ܺ�CR=0.5��lena.png','png'); 
load('Encryption_Key.mat');

U1=I(:,:,1);  
U2=I(:,:,2);  
U3=I(:,:,3);  

[M,N]=size(U1);                      %��ͼ������и�ֵ��M,N

%��Կ���жȷ���
t=KEY(1,1);    %�ֿ��С

SUM=M*N;

%��Կ���жȷ���
% u=3.9999000000000001;
u=KEY(1,2);
% u=3.78774868768;

xx0=KEY(1,10);
% xx0=0.3884000000000001;

xx1=KEY(1,11);
% xx1=0.4133000000000001;

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
% %% 2.����Logistic��������
% % u=3.990000000000001; %��Կ�����Բ���  10^-15
% % u=3.99;%��Կ��Logistic������
% % x0=0.7067000000000001; %��Կ�����Բ���  10^-16
% % x0=0.5475; %��Կ��Logistic��ֵx0
%  
% x0=KEY(1,3);
% % x0=0.547623608200001;
% % x0=0.8477234627836;
% 
% p=zeros(1,SUM+1000);
% p(1)=x0;
% for i=1:SUM+999                        %����SUM+999��ѭ��������SUM+1000������
%     p(i+1)=u*p(i)*(1-p(i));
% end
% p=p(1001:length(p));
% 
% %% 3.��p���б任��0~255��Χ��������ת����M*N�Ķ�ά����R
% p=mod(round(p*10^4),256);
% R=reshape(p,N,M)';  %ת��M��N��
% 
% %% 4.�����緽��
% % ���ĸ���ֵX0,Y0,Z0,H0
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
% %%������Chenϵͳ
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
% %%�����׳�����Chenϵͳ
% % h=0.01;%���ò���
% % NN=25000;%���л���ϵͳ�õ������и���
% % q=0.95;%�����׽���
% % z0=[X0 Y0 Z0 H0];
% % [~,y]=FrataSim(h,NN,z0,q);
% % X=y(1,:);
% % X=X(5000:5000+r);        %ȥ��ǰ3001���ø��õ�����ԣ�������ϵͳ���Ӻ����������3000�㣩
% % Y=y(2,:);
% % Y=Y(5000:5000+r);
% % Z=y(3,:);
% % Z=Z(5000:5000+r);
% % H=y(4,:);
% % % H=H(5000:length(H));
% % H=H(5000:5000+r);
% 
% %��Գ���512*512�ĸ����ͼ����Ҫ����Ļ�������
% h=0.01;%���ò���
% NN=35000;%���л���ϵͳ�õ������и���
% q=0.95;%�����׽���
% 
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
% %% 5.DNA����
% %X,Y�ֱ����I��R��DNA���뷽ʽ����8�֣�1~8
% X=mod(round(X*10^4),8)+1;
% Y=mod(round(Y*10^4),8)+1;
% Z=mod(round(Z*10^4),4);
% Z(Z==0)=4;      %�Ӽ�������
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
%     Q2_Y=DNA_yunsuan(Q1_Y,Q1_last_Y,Z(i));        %��ɢǰ
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
% %% �����ݴ���
% rate=4;     %ѹ������1���ڵ�������������Ϊ��Կ
% max_Multiple_first=KEY(1,12);%������Կ֮һ
% % max_Multiple_first=9999;%������Կ֮һ
% 
% max_Multiple=KEY(1,13);%������Կ֮һ
% % max_Multiple=666;
% 
% load('Mask.mat');
% load('Coefficient_Matrix.mat');
% T1=zeros(t,t);  % ѹ������
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
% % T1=T1.*YS1;%û��ѹ�������
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
% fun = @(block_struct) block_struct.data./T1;         %Matlab������blockproc��������blkproc����
% J1 = blockproc(U1,[t t],fun);                        %�������ݴ���֮��ľ���
% J2 = blockproc(U2,[t t],fun);                        %�������ݴ���֮��ľ���
% J3 = blockproc(U3,[t t],fun);                        %�������ݴ���֮��ľ���
% %% ��DCT����
% T=dctmtx(t);        %8��DCT�任����
% fun = @(block_struct) T'*block_struct.data*T;       %Matlab������blockproc��������blkproc����
% K1 = blockproc(J1,[t t],fun);              %DCT��任���ع�ͼ��
% K2 = blockproc(J2,[t t],fun);              %DCT��任���ع�ͼ��
% K3 = blockproc(J3,[t t],fun);              %DCT��任���ع�ͼ��
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
% %% 6��ȥ������ʱ������
% M1=KEY(1,8);   %����ʱ����Ĳ�����M1=mod(M,t);��Ϊ��Կ
% N1=KEY(1,9);   %����ʱ����Ĳ�����N1=mod(N,t);��Ϊ��Կ
% if M1~=0
%     I_jiemi=I_jiemi(1:M-t+M1,:,:);
% end
% if N1~=0
%     I_jiemi=I_jiemi(:,1:N-t+N1,:);
% end
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,1));
% %title('����ͼƬRͨ��ֱ��ͼ');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/lena/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.1.05����/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬRͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/����ͼƬRͨ��ֱ��ͼ.png','png');
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,2));
% % title('����ͼƬGͨ��ֱ��ͼ');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/lena/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.1.05����/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬGͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/����ͼƬGͨ��ֱ��ͼ.png','png');
% 
% B=figure('color',[1 1 1]);
% imhist(I_jiemi(:,:,3));
% % title('����ͼƬBͨ��ֱ��ͼ');
% frame = getframe(B);
% B_jiamizhifangtu=frame2im(frame);
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/lena/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.1.05����/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/����/����ͼƬBͨ��ֱ��ͼ.png','png');
% % imwrite(B_jiamizhifangtu,'../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/����ͼƬBͨ��ֱ��ͼ.png','png');
% 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/���ܺ��lena.png','png'); 
% 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/u+10^-16.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/X(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/Y(0)+10^-16.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/Z(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/H(0)+10^-14.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/x0+10^-15.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/x1+10^-16.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/x2+10^-16.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/p+10^-16.png','png'); 
% 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/lena/all_dontknow.png','png');
% 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��peppers.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��baboon.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/�ɻ�/���ܺ��airplane.png','png'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/4.1.05����/���ܺ��4.1.05.tiff','tiff'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/���ܺ��4.2.06.tiff','tiff'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/����/���ܺ��house.tiff','tiff'); 
% % imwrite(I_jiemi,'../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/���ܺ��kodim08.png','png'); 
% 
% disp('������Ľ�����ԿΪ��');
% disp(['��Կ1��t=',num2str(t),'     ��Կ2��u=',num2str(u),'    ��Կ3��x0=',num2str(x0),'    ��Կ4��x(0)=',num2str(X0),'   ��Կ5��y(0)=',num2str(Y0) '   ��Կ6��z(0)=',num2str(Z0)]);
% disp(['��Կ7��h(0)=',num2str(H0),'   ��Կ8��M1=',num2str(M1),'   ��Կ9��N1=',num2str(N1),'   ��Կ10��xx0=',num2str(xx0),'   ��Կ11��xx1=',num2str(xx1) '   ��Կ12��max_Multiple_first=',num2str(max_Multiple_first) '   ��Կ13��max_Multiple=',num2str(max_Multiple)]);
% disp('�������'); 
% figure;imshow(I_jiemi);
% title('���ܺ�ͼƬ');
% % figure;imshow(V1);
% % title('���ܺ�ͼƬ2');