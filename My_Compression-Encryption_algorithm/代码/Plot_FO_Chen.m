%%�������׳�����Chenϵͳ��������ͼ
%%ʱ�䣺5.11.2020
%%���ߣ����

clear;clc;
I=imread('../ԭʼ�����ܡ�����ͼƬ/lena/lena.png','png');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/peppers.png','png');       %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/baboon.png','png');         %��ȡͼ����Ϣ������
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�ɻ�/airplane.png','png');         %��ȡͼ����Ϣ���ɻ�
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.1.05����/4.1.05.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/4.2.06С��/4.2.06.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/����/house.tiff','tiff');         %��ȡͼ����Ϣ
% I=imread('../ԭʼ�����ܡ�����ͼƬ/�Ǳ�/kodim08.png','png');         %��ȡͼ����Ϣ
figure;imshow(I);title('ԭʼͼƬ');

I1=I(:,:,1);        %R
I2=I(:,:,2);        %G
I3=I(:,:,3);        %B
[M,N]=size(I1);                      %��ͼ������и�ֵ��M,N
t=4;    %�ֿ��С
%% 1.����
%��ͼ��������������ɿ��Ա�t����������tΪ�ֿ�Ĵ�С��
M1=mod(M,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0
N1=mod(N,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0
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
%% 4.��������Chen�ϳ�����ϵͳ
%���ĸ���ֵX0,Y0,Z0,H0
r=(M/t)*(N/t);      %rΪ�ֿ����
%����ĸ���ֵ
X0=sum(sum(bitand(I1,17)))/(17*SUM);
Y0=sum(sum(bitand(I2,34)))/(34*SUM);
Z0=sum(sum(bitand(I3,68)))/(68*SUM);
H0=sum(sum(bitand(I1,136)))/(136*SUM);
%������λС��
X0=round(X0*10^4)/10^4;%��Ϊ��Կ3
Y0=round(Y0*10^4)/10^4;%��Ϊ��Կ4
Z0=round(Z0*10^4)/10^4;%��Ϊ��Կ5
H0=round(H0*10^4)/10^4;%��Ϊ��Կ6
%���ݳ�ֵ�����Chen�ϳ�����ϵͳ���õ��ĸ���������  
%������Chenϵͳ
% A=chen_output(X0,Y0,Z0,H0,r);
% X=A(:,1);
% X=X(3002:length(X));
% Y=A(:,2);
% Y=Y(3002:length(Y));
% Z=A(:,3);
% Z=Z(3002:length(Z));
% H=A(:,4);
% H=H(3002:length(H));

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

%%����������ͼ
Q=figure(1);
plot(X,Y);
xlabel('\itx')
ylabel('\ity')
frame = getframe(Q);
x_y=frame2im(frame);
imwrite(x_y,'../����ͼ/lena/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/����/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/����/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/�ɻ�/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/4.1.05����/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/4.2.06С��/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/����/�����׳�����Chenϵͳ������x-y.png','png');
% imwrite(x_y,'../����ͼ/�Ǳ�/�����׳�����Chenϵͳ������x-y.png','png');

W=figure(2);
plot(X,Z)
xlabel('\itx')
ylabel('\itz')
frame = getframe(W);
x_z=frame2im(frame);
imwrite(x_z,'../����ͼ/lena/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/����/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/����/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/�ɻ�/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/4.1.05����/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/4.2.06С��/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/����/�����׳�����Chenϵͳ������x-z.png','png');
% imwrite(x_z,'../����ͼ/�Ǳ�/�����׳�����Chenϵͳ������x-z.png','png');

E=figure(3);
plot(X,H)
xlabel('\itx')
ylabel('\itw')
frame = getframe(E);
x_w=frame2im(frame);
imwrite(x_w,'../����ͼ/lena/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/����/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/����/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/�ɻ�/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/4.1.05����/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/4.2.06С��/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/����/�����׳�����Chenϵͳ������x-w.png','png');
% imwrite(x_w,'../����ͼ/�Ǳ�/�����׳�����Chenϵͳ������x-w.png','png');

R=figure(4);
plot(Z,H)
xlabel('\itz')
ylabel('\itw')
frame = getframe(R);
z_w=frame2im(frame);
imwrite(z_w,'../����ͼ/lena/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/����/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/����/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/�ɻ�/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/4.1.05����/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/4.2.06С��/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/����/�����׳�����Chenϵͳ������z-w.png','png');
% imwrite(z_w,'../����ͼ/�Ǳ�/�����׳�����Chenϵͳ������z-w.png','png');