clear ;
clc;
close all;
%% ��������ݣ�����1.��ͼA,b 2.��ͨ���� mu ,Sigma 3.�û����� zeta 4.�㷨��ֹ����num��gap

start_node=[1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,8,8,8,8,9,9,9,10,10,10,10,...
    10,11,11,11,11,12,12,12,13,13,14,14,14,15,15,15,15,16,16,16,16,17,17,...
    17,18,18,18,19,19,19,20,20,20,20,21,21,21,22,22,22,22,23,23,23,24,24,24];
end_node=[2,3,1,6,1,4,12,3,5,11,4,6,9,2,5,8,8,18,6,7,9,16,5,8,10,9,11,15,...
    16,17,4,10,12,14,3,11,13,12,24,11,15,23,10,14,19,22,8,10,17,18,10,16,...
    19,7,16,20,15,17,20,18,19,21,22,20,22,24,15,20,21,23,14,22,24,13,21,23];
A=zeros(24,76);
for i=1:76
    A(start_node(i),i)=1;
end
for i=1:76
    A(end_node(i),i)=-1;
end
mu=10*randi([1,8],1,76)';

for i=1:76
    standard_deviation(i)=6^3*5/mu(i);
end
covSigma=zeros(76);
for i=1:76
    covSigma(i,i)=standard_deviation(i)^2;
end
b=zeros(24,1)';
b(1)=1;

number_zeta=100;
zeta=zeros(number_zeta,1);
for i=1:number_zeta
    zeta(i)=i;
end% zeta is tolerance
%% �յ��node2��node24
number_end_node=23;

b=zeros(24,1);
b(1)=1;
b(15)=-1;
gap=0.01;
%% zeta ��1��100

ZETA=zeta(6);

[path_rsp_MIQP]=func_mixed_integer_quadradic_programming(A,b,covSigma,mu,gap,ZETA);


