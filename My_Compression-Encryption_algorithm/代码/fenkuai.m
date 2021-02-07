%% 子函数 分块函数，t为每块的边长，I为要分块的图像，num为返回第几大块,数据块从左向右依次取下来，然后再往下取
function fv=fenkuai(t,I,num)
[~,N]=size(I);
S=N/t;
x=floor(num/S)+1;      %第几大行
y=mod(num,S);           %第几大列
if y==0
    x=x-1;
    y=S;
end
fv=I(t*(x-1)+1:t*x,t*(y-1)+1:t*y);


