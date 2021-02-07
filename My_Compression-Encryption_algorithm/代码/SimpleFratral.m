%%子函数，求解分数阶超混沌Chen系统
%利用Adomian分解方法进行求解
%时间：2020.5.15
%作者：董昊
function dy = SimpleFratral(z0,h,q)
dy = zeros(1,4);
if q<=1&&q>0
    c10=z0(1,1);
    c20=z0(1,2);
    c30=z0(1,3);
    c40=z0(1,4);
end
if q>1&&q<0
    c10=z0(1)+z0(1)*h;
    c20=z0(2)+z0(2)*h;
    c30=z0(3)+z0(3)*h;
    c40=z0(4)+z0(4)*h;
end

c11=35*(c20-c10)+c40;
c21=7*c10-c10*c30+12*c20;
c31=c10*c20-3*c30;
c41=c20*c30+0.5*c40;

c12=35*(c21-c11)+c41;
c22=7*c11+12*c21-c11*c30-c10*c31;
c32=-3*c31+c11*c20+c10*c21;
c42=0.5*c41+c21*c30+c20*c31;

c13=35*(c22-c12)+c42;
c23=7*c12+12*c22-c12*c30-c10*c32-c11*c31*(gamma(2*q+1)/gamma(q+1)^2);
c33=-3*c32+c12*c20+c10*c22+c11*c21*(gamma(2*q+1)/gamma(q+1)^2);
c43=0.5*c42+c22*c30+c20*c32+c21*c31*(gamma(2*q+1)/gamma(q+1)^2);

c14=35*(c23-c13)+c43;
c24=7*c13+12*c23-c13*c30-c10*c33-(c11*c32+c12*c31)*(gamma(3*q+1)/(gamma(2*q+1)*gamma(q+1)));
c34=-3*c33+c13*c20+c10*c23+(c11*c22+c12*c21)*(gamma(3*q+1)/(gamma(2*q+1)*gamma(q+1)));
c44=0.5*c43+c23*c30+c20*c33+(c22*c31+c21*c32)*(gamma(3*q+1)/(gamma(2*q+1)*gamma(q+1)));

c15=35*(c24-c14)+c44;
c25=7*c14+12*c24-c14*c30-c10*c34-(c13*c31+c11*c33)*(gamma(4*q+1)/(gamma(3*q+1)*gamma(q+1)))-c12*c32*(gamma(4*q+1)/gamma(2*q+1)^2);
c35=-3*c34+c14*c20+c10*c24+(c13*c21+c11*c23)*(gamma(4*q+1)/(gamma(3*q+1)*gamma(q+1)))+c12*c22*(gamma(4*q+1)/gamma(2*q+1)^2);
c45=0.5*c44+c24*c30+c20*c34+(c23*c31+c21*c33)*(gamma(4*q+1)/(gamma(3*q+1)*gamma(q+1)))+c22*c32*(gamma(4*q+1)/gamma(2*q+1)^2);

c16=35*(c25-c15)+c45;
c26=7*c15+12*c25-c15*c30-c10*c35-(c14*c31+c11*c34)*(gamma(5*q+1)/(gamma(4*q+1)*gamma(q+1)))-(c13*c32+c12*c33)*(gamma(5*q+1)/(gamma(3*q+1)*gamma(2*q+1)));
c36=-3*c35+c15*c20+c10*c25+(c14*c21+c11*c24)*(gamma(5*q+1)/(gamma(4*q+1)*gamma(q+1)))+(c13*c22+c12*c23)*(gamma(5*q+1)/(gamma(3*q+1)*gamma(2*q+1)));
c46=0.5*c45+c25*c30+c20*c35+(c24*c31+c21*c34)*(gamma(5*q+1)/(gamma(4*q+1)*gamma(q+1)))+(c23*c32+c22*c33)*(gamma(5*q+1)/(gamma(3*q+1)*gamma(2*q+1)));

dy(1)=c10+c11*(h^q/gamma(q+1))+c12*h^(2*q)/gamma(2*q+1)+c13*h^(3*q)/gamma(3*q+1)+c14*h^(4*q)/gamma(4*q+1)+c15*h^(5*q)/gamma(5*q+1)+c16*h^(6*q)/gamma(6*q+1);
dy(2)=c20+c21*(h^q/gamma(q+1))+c22*h^(2*q)/gamma(2*q+1)+c23*h^(3*q)/gamma(3*q+1)+c24*h^(4*q)/gamma(4*q+1)+c25*h^(5*q)/gamma(5*q+1)+c26*h^(6*q)/gamma(6*q+1);
dy(3)=c30+c31*(h^q/gamma(q+1))+c32*h^(2*q)/gamma(2*q+1)+c33*h^(3*q)/gamma(3*q+1)+c34*h^(4*q)/gamma(4*q+1)+c35*h^(5*q)/gamma(5*q+1)+c36*h^(6*q)/gamma(6*q+1);
dy(4)=c40+c41*(h^q/gamma(q+1))+c42*h^(2*q)/gamma(2*q+1)+c43*h^(3*q)/gamma(3*q+1)+c44*h^(4*q)/gamma(4*q+1)+c45*h^(5*q)/gamma(5*q+1)+c46*h^(6*q)/gamma(6*q+1);

























% function fy=FO_Chen_output(x0,y0,z0,h0)
% a=35;
% b=7;
% c=12;
% d=3;
% r=0.5;
% fdefun = @(t,y)[a*(y(2)-y(1))+y(4);
%                 b*y(1)-y(1)*y(3)+c*y(2);
%                 y(1)*y(2)-d*y(3);
%                 y(2)*y(3)+r*y(4)];
% alpha = 0.95;
% t0=0;
% tfinal=500;
% y0=[x0;y0;z0;h0] ;
% h=2^(-7);
% % h=0.01;
% [~, y_fde12] = fde12(alpha,fdefun,t0,tfinal,y0,h);
% fy=y_fde12;
% end

%%做混沌系统的吸引子图
% figure(1)
% plot(t,y_fde12(1,:),t,y_fde12(2,:),t,y_fde12(3,:),t,y_fde12(4,:));
% xlabel('t');ylabel('y(t)');
% legend('y_1(t)','y_2(t)','y_3(t)','y_4(t)');
% title('FDE solved by the efd12.m code');
% figure
% plot(y_fde12(1,:),y_fde12(2,:))
% xlabel('\itx')
% ylabel('\ity')
% figure
% plot(y_fde12(1,:),y_fde12(3,:))
% xlabel('\itx')
% ylabel('\itz')
% figure
% plot(y_fde12(1,:),y_fde12(4,:))
% xlabel('\itx')
% ylabel('\itw')
% figure
% plot(y_fde12(3,:),y_fde12(4,:))
% xlabel('\itz')
% ylabel('\itw')
% zlabel('\itz')

% %%测试求解分数阶微分方程组
% %example
% a = 1; mu = 4;
% fdefun = @(t,y)[a-(mu+1)*y(1)+y(1)^2*y(2) ; mu*y(1)-y(1)^2*y(2)];
% alpha = 0.8 ;
% t0 = 0 ; tfinal = 100 ; y0 = [ 0.2 ; 0.03] ;
% h = 2^(-6) ;
% [t, y_fde12] = fde12(alpha,fdefun,t0,tfinal,y0,h);
% figure(1)
% plot(t,y_fde12(1,:),t,y_fde12(2,:));
% xlabel('t');ylabel('y(t)');
% legend('y_1(t)','y_2(t)');
% title('FDE solved by the efd12.m code');