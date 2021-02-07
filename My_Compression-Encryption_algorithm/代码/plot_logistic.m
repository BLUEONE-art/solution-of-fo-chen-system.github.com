clc;  
clear all; 
close all;

u=0:0.001:4;
a=length(u);
x=0.4*ones(1,4001);
N1=400;
N2=100;
f=zeros(N1+N2,length(u));
for i = 1:N1+N2
   x=u.*x.*(1-x); 
   f(i,:)=x;
end
f=f(N1+1:end,:);
plot(u,f,'k.','MarkerSize',0.5);   
title('\fontsize{10}Logisticӳ�����ͼ');          
xlabel('\fontsize{10}��֧����u'),
ylabel('\fontsize{10}������зֲ�x(n)'); 
