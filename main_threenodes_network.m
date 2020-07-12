clear;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap
A=zeros(3,40);
A(1,1:20)=1;
A(2,1:20)=-1;
A(2,21:40)=1;
A(3,21:40)=-1;
%A = generate_graph();
b=[1,0,-1]';    % b is the OD column vector

%% generate mean value and variance of edges
mu=zeros(40,1);
for i=1:20
    mu(i)=i;
end
for i=21:40
    mu(i)=i-20;
end
covSigma=zeros(40,40);
for i=1:40
    covSigma(i,i)=(100/mu(i))^2;
end
%% generate zeta from 1 to 100
number_zeta=100;
zeta=zeros(number_zeta,1);
for i=1:number_zeta
    zeta(i)=i;
end% zeta is tolerance
gap = 0.01;   % set the value of gap
max_iter=100;
performance_total=zeros(number_zeta,400);
performace_optimal_value=zeros(number_zeta,1);
performance_GP3=zeros(number_zeta,1);
counter=1;
path_name=zeros(40,400);
for j=1:20
    for z=1:20
        path_name(j,counter)=1;
        path_name(z+20,counter)=1;
        counter=counter+1;
    end
end

true_gap = zeros(number_zeta,1);
for i=1:1%% only test the case: zeta = 18

    for j=1:400
        performance_total(i,j)=mu'*path_name(:,j)+zeta(i)*sqrt(path_name(:,j)'*covSigma*path_name(:,j));
    end
    performace_optimal_value(i)=min(performance_total(i,:));
    [path_rsp_GP3,iter_GP3,gap_GP3(i),UB_GP3,LB_GP3,gap_total_GP3,gap_total_real_GP3,iteration_GP3]=func_GP3(A,b,mu,covSigma,zeta(i),gap);

    performance_GP3(i)=mu'*path_rsp_GP3+zeta(i)*sqrt(path_rsp_GP3'*covSigma*path_rsp_GP3);
    true_gap(i)= performance_GP3(i)-performace_optimal_value(i);
end
figure(1)
plot(zeta,performace_optimal_value,zeta,performance_GP3);
legend('optimal','GP3');
xlabel('zeta');
ylabel('the value of performance');
title('performance');


