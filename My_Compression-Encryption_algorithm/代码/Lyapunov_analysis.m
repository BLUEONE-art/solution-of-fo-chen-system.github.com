%%计算混沌系统吸引子的Lyapunov指数
%%作者：董昊
%%时间：2020.05.11
clear;
yinit=[1,1,1,1];
orthmatrix=[1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 1];
a=35;
b=7;
c=12;
d=3;
alpha = 0.95;
r=0.5;
% fdefun = @(t,y)[a*(y(2)-y(1))+y(4);
%                 b*y(1)-y(1)*y(3)+c*y(2);
%                 y(1)*y(2)-d*y(3);
%                 y(2)*y(3)+r*y(4)];
y=zeros(20,1);
%初始化输入
y(1:4)=yinit;
y(5:20)=orthmatrix;
tstart=0;%时间初始值
tstep=1e-3;%时间步长
wholetimes=1e5;%总的循环次数
steps=10;%每次演化的步数
iteratetimes=wholetimes/steps;%演化的次数?
mod=zeros(4,1);
lp=zeros(4,1);
%初始化四个Lyapunov指数
Lyapunov1=zeros(iteratetimes,1);
Lyapunov2=zeros(iteratetimes,1);
Lyapunov3=zeros(iteratetimes,1);
Lyapunov4=zeros(iteratetimes,1);
for i=1:iteratetimes
    tspan=tstart:tstep:(tstart+tstep*steps);
    [T,Y]=ode45('fra_chaos_fun',tspan,y);
%     [T1,Y1] = fde12(alpha,fdefun,tstart,(tstart+tstep*steps),h(1:4),tstep);
%     [T2,Y2] = fde12(alpha,fdefun,tstart,(tstart+tstep*steps),h(5:8),tstep);
%     [T3,Y3] = fde12(alpha,fdefun,tstart,(tstart+tstep*steps),h(9:12),tstep);
%     [T4,Y4] = fde12(alpha,fdefun,tstart,(tstart+tstep*steps),h(13:16),tstep);
%     [T5,Y5] = fde12(alpha,fdefun,tstart,(tstart+tstep*steps),h(17:20),tstep);
%     Y=[Y1;Y2;Y3;Y4;Y5];
    %取积分得到的最后一个时刻的值
%     Y=Y';
    y=Y(size(Y,1),:);
    %重新定义起始时刻
    tstart=tstart + tstep*steps;
    y0=[y(5) y(9)  y(13) y(17);
        y(6) y(10) y(14) y(18);
        y(7) y(11) y(15) y(19);
        y(8) y(12) y(16) y(20)];
    %正交化
    y0=ThreeGS(y0);
    %取四个向量的模
    mod(1)=sqrt(y0(:,1)'*y0(:,1));
    mod(2)=sqrt(y0(:,2)'*y0(:,2));
    mod(3)=sqrt(y0(:,3)'*y0(:,3));
    mod(4)=sqrt(y0(:,4)'*y0(:,4));
    y0(:,1)=y0(:,1)/mod(1);
    y0(:,2)=y0(:,2)/mod(2);
    y0(:,3)=y0(:,3)/mod(3);
    y0(:,4)=y0(:,4)/mod(4);
    lp=lp+log(abs(mod));
    %三个Lyapunov指数
    Lyapunov1(i)=lp(1)/(tstart);
    Lyapunov2(i)=lp(2)/(tstart);
    Lyapunov3(i)=lp(3)/(tstart);
    Lyapunov4(i)=lp(4)/(tstart);
    y(5:20)=y0';
    y=y';
end
%作Lyapunov指数谱图
i=1:iteratetimes;
plot(i,Lyapunov1)
hold on
plot(i,Lyapunov2)
hold on
plot(i,Lyapunov3)
hold on
plot(i,Lyapunov4)
xlabel('迭代次数');ylabel('Lyapunov指数');
legend('Lyapunov1','Lyapunov2','Lyapunov3','Lyapunov4');
title('分数阶超混沌Chen系统Lyapunov指数');

% %%计算混沌系统吸引子的Lyapunov指数
% %%作者：董昊
% %%时间：2020.05.11
% clear;
% yinit=[1,1,1];
% orthmatrix=[1 0 0;
%             0 1 0;
%             0 0 1];
% a=0.15;
% b=0.20;
% c=10.0;
% y=zeros(12,1);
% %初始化输入
% y(1:3)=yinit;
% y(4:12)=orthmatrix;
% tstart=0;%时间初始值
% tstep=1e-3;%时间步长
% wholetimes=1e5;%总的循环次数
% steps=10;%每次演化的步数
% iteratetimes=wholetimes/steps;%演化的次数?
% mod=zeros(3,1);
% lp=zeros(3,1);
% %初始化三个Lyapunov指数
% Lyapunov1=zeros(iteratetimes,1);
% Lyapunov2=zeros(iteratetimes,1);
% Lyapunov3=zeros(iteratetimes,1);
% for i=1:iteratetimes
%     tspan=tstart:tstep:(tstart+tstep*steps);
%     [T,Y]=ode45('Rossler_ly',tspan,y);
%     %取积分得到的最后一个时刻的值
%     y=Y(size(Y,1),:);
%     %重新定义起始时刻
%     tstart=tstart + tstep*steps;
%     y0=[y(4) y(7) y(10);
%         y(5) y(8) y(11);
%         y(6) y(9) y(12)];
%     %正交化
%     y0=ThreeGS(y0);
%     %取三个向量的模
%     mod(1)=sqrt(y0(:,1)'*y0(:,1));
%     mod(2)=sqrt(y0(:,2)'*y0(:,2));
%     mod(3)=sqrt(y0(:,3)'*y0(:,3));
%     y0(:,1)=y0(:,1)/mod(1);
%     y0(:,2)=y0(:,2)/mod(2);
%     y0(:,3)=y0(:,3)/mod(3);
%     lp=lp+log(abs(mod));
%     %三个Lyapunov指数
%     Lyapunov1(i)=lp(1)/(tstart);
%     Lyapunov2(i)=lp(2)/(tstart);
%     Lyapunov3(i)=lp(3)/(tstart);
%     y(4:12)=y0';
% end
% %作Lyapunov指数谱图
% i=1:iteratetimes;
% plot(i,Lyapunov1,Lyapunov2,Lyapunov3)








