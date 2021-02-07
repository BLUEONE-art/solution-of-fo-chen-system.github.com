clc;clear;
h=0.01;
NN=10000;
q=0.95;%分数阶阶数
z0=[0.1 0.2 0.3 0.4];%初值
[t,y]=FrataSim(h,NN,z0,q);
R=figure('color',[1 1 1]);
plot(y(1,:),y(2,:));
xlabel('\itx');
ylabel('\ity');
frame = getframe(R);
x_y=frame2im(frame);
imwrite(x_y,'../混沌图/测试/分数阶超混沌Chen系统吸引子x-y.png','png');

R=figure('color',[1 1 1]);
plot(y(1,:),y(3,:));
xlabel('\itx');
ylabel('\itz');
frame = getframe(R);
x_z=frame2im(frame);
imwrite(x_z,'../混沌图/测试/分数阶超混沌Chen系统吸引子x_z.png','png');

R=figure('color',[1 1 1]);
plot(y(1,:),y(4,:));
xlabel('\itx');
ylabel('\itw');
frame = getframe(R);
x_w=frame2im(frame);
imwrite(x_w,'../混沌图/测试/分数阶超混沌Chen系统吸引子x_w.png','png');

R=figure('color',[1 1 1]);
plot(y(3,:),y(4,:));
xlabel('\itz');
ylabel('\itw');
frame = getframe(R);
z_w=frame2im(frame);
imwrite(z_w,'../混沌图/测试/分数阶超混沌Chen系统吸引子z_w.png','png');

% %%G-S正交化
% function A=ThreeGS(V)%为3*3的向量
% v1=V(:,1);
% v2=V(:,2);
% v3=V(:,3);
% v4=V(:,4);
% a1=zeros(4,1);
% a2=zeros(4,1);
% a3=zeros(4,1);
% a4=zeros(4,1);
% a1=v1;
% a2=v2-((a1'*v2)/(a1'*a1))*a1;
% a3=v3-((a1'*v3)/(a1'*a1))*a1-((a2'*v3)/(a2'*a2))*a2;
% a4=v4-((a1'*v4)/(a1'*a1))*a1-((a2'*v4)/(a2'*a2))*a2-((a3'*v4)/(a3'*a3))*a3;
% A=[a1,a2,a3,a4];
% % A=[a1,a2,a3];

