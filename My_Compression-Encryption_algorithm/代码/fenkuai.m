%% �Ӻ��� �ֿ麯����tΪÿ��ı߳���IΪҪ�ֿ��ͼ��numΪ���صڼ����,���ݿ������������ȡ������Ȼ��������ȡ
function fv=fenkuai(t,I,num)
[~,N]=size(I);
S=N/t;
x=floor(num/S)+1;      %�ڼ�����
y=mod(num,S);           %�ڼ�����
if y==0
    x=x-1;
    y=S;
end
fv=I(t*(x-1)+1:t*x,t*(y-1)+1:t*y);


