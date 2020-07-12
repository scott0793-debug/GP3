function [ mu,Sigma,RHO ] = covarianceMatrix( inputData )
% ? covarianceMatrix( inputData )
% ? ����һ������Э�������ĺ���
% ? inputData ? ��������
% ? ÿһ��Ϊһ��ά��
% ? ÿһ��Ϊһ������
%�����������ά��
[m,n] = size(inputData);
%����Э�������
Sigma = zeros(m,m);
%ȡ��ÿά����ƽ��ֵ
mu = zeros(m,1);
for i = 1:m
    mu(i) = mean(inputData(i,:));
end
%����Э����
for i = 1:m
    for j = 1:m
        Sigma(i,j) = ((inputData(i,:)-mu(i))*(inputData(j,:)-mu(j))')./(n-1);
    end
end
%�������ϵ��
RHO = corr(inputData');

