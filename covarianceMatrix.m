function [ mu,Sigma,RHO ] = covarianceMatrix( inputData )
% ? covarianceMatrix( inputData )
% ? 这是一个计算协方差矩阵的函数
% ? inputData ? 输入数据
% ? 每一行为一个维度
% ? 每一列为一个样本
%获得输入数据维度
[m,n] = size(inputData);
%创建协方差矩阵
Sigma = zeros(m,m);
%取得每维数据平均值
mu = zeros(m,1);
for i = 1:m
    mu(i) = mean(inputData(i,:));
end
%计算协方差
for i = 1:m
    for j = 1:m
        Sigma(i,j) = ((inputData(i,:)-mu(i))*(inputData(j,:)-mu(j))')./(n-1);
    end
end
%计算相关系数
RHO = corr(inputData');

