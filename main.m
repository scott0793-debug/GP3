
clear ;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap
A=[1,-1,0;1,0,-1;0,1,-1;0,1,-1]'; % A is the map matrix
%A = generate_graph();
b=[1,0,-1]';    % b is the OD column vector

inputData=[15,10,5;21,20,19;2.2,10.1,18;18,10.1,2.2];% 每条路径的行驶时间数据
% generate mean  value E, covariance matrix COVMAT, correlation coefficient maxtrix RHO
[ mu,Sigma,RHO ] = covarianceMatrix( inputData );

%% 将协方差矩阵改为正定矩阵covSigma
err_cnt = 0;
for i = 1:1000
    try 
        [m,n]=size(Sigma);
        covSigma = cov(Sigma) + .0001 * eye(m); 
        m = mean(Sigma); 
        mvnpdf(Sigma, m, covSigma);
    catch me
        err_cnt = err_cnt + 1;
    end
end
clear err_ent

zeta = 0.4;      % zeta is tolerance
gap = 0.01;   
max_iter=100;

[x]=func_GP3(A,b,mu,covSigma,zeta,gap);
path_rsp_zwl=func_rsp_zwl(A,b,mu,covSigma,zeta,max_iter,gap);
path_rsp_zyf=func_rsp_zyf(A,b,mu,covSigma,zeta,max_iter,gap);
path_rsp_CG=func_rsp_CG(A,b,mu,covSigma,zeta,max_iter,gap);
