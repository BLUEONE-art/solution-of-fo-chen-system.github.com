%%子函数，求解分数阶超混沌Chen系统
%利用Adomian分解方法进行求解
%时间：2020.5.15
%作者：董昊
function [t,y]=FrataSim(h,NN,z0,q)
y=zeros(4,NN);
t=zeros(1,NN);
tt=0;
for i=1:NN
    t(i)=tt;
    dy=SimpleFratral(z0,h,q);
    y(:,i)=dy;
    z0=dy;
    tt=tt+h;
end