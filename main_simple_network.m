clear ;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap
A=zeros(2,50);
A(1,:)=1;
A(2,:)=-1;
%A = generate_graph();
b=[1,-1]';    % b is the OD column vector

mu=zeros(50,1);
for i=1:50
    mu(i)=i;
end
covSigma=zeros(50);
for i=1:50
    covSigma(i,i)=(100/i)^2;
end

number_zeta=100;
zeta=zeros(number_zeta,1);
for i=1:number_zeta
    zeta(i)=i;
end% zeta is tolerance
gap = 0.01;   
max_iter=100;
performance_total=zeros(number_zeta,50);
%for i=1:number_zeta
for i=1:100
    for j=1:50
        path_name=zeros(50,1);
        path_name(j)=1;
        performance_total(i,j)=mu'*path_name+zeta(i)*sqrt(path_name'*covSigma*path_name);
    end
    performace_optimal_value(i)=min(performance_total(i,:));
    tic;
    [path_rsp_zyf,iter_zyf,gap_zyf(i),upper_bound_total_zyf,lower_bound_total_zyf,gap_total_zyf,gap_real_total_zyf,iteration_zyf]=func_rsp_zyf(A,b,mu,covSigma,zeta(i),max_iter,gap);
    t_zyf(i)=toc;
    performance_zyf(i)=mu'*path_rsp_zyf+zeta(i)*sqrt(path_rsp_zyf'*covSigma*path_rsp_zyf);
    tic;
    [path_rsp_zwl,iter_zwl,gap_zwl(i),UB_zwl,LB_zwl,gap_total_zwf,gap_LB_UB_real,iteration_zwl]=func_rsp_zwl(A,b,mu,covSigma,zeta(i),max_iter,gap);
    t_zwl(i)=toc;
    performance_zwl(i)=mu'*path_rsp_zwl+zeta(i)*sqrt(path_rsp_zwl'*covSigma*path_rsp_zwl);
    tic;
    [path_rsp_SP,iter_SP,gap_RDG(i),upper_bound_total_SP,lower_bound_total_SP,gap_RDG_total_SP,gap_RDG_real_SP,iteration_SP]=func_rsp_SP(A,b,mu,covSigma,zeta(i),max_iter,gap);
    t_SP(i)=toc;
    performance_SP(i)=mu'*path_rsp_SP+zeta(i)*sqrt(path_rsp_SP'*covSigma*path_rsp_SP);
    tic;
    [path_rsp_GP3,iter_GP3,gap_GP3(i),UB_GP3,LB_GP3,gap_total_GP3,gap_total_real_GP3,iteration_GP3]=func_GP3(A,b,mu,covSigma,zeta(i),gap);
    t_GP3(i)=toc;
    performance_GP3(i)=mu'*path_rsp_GP3+zeta(i)*sqrt(path_rsp_GP3'*covSigma*path_rsp_GP3);
    table_GP3=table(iteration_GP3',UB_GP3',LB_GP3',gap_total_GP3',gap_total_real_GP3','VariableNames',{'iteration','UB','LB','gap','real_gap'});
    table_SP=table(iteration_SP',upper_bound_total_SP',lower_bound_total_SP',gap_RDG_total_SP',gap_RDG_real_SP','VariableNames',{'iteration','UB','LB','gap','real_gap'});
    table_zwl=table(iteration_zwl',UB_zwl',LB_zwl',gap_total_zwf',gap_LB_UB_real','VariableNames',{'iteration','UB','LB','gap','real_gap'});
    table_zyf=table(iteration_zyf',upper_bound_total_zyf',lower_bound_total_zyf',gap_total_zyf',gap_real_total_zyf','VariableNames',{'iteration','UB','LB','gap','real_gap'});
end
figure(1)
plot(zeta,gap_GP3,zeta,gap_zyf,zeta,gap_RDG,zeta,gap_zwl);
legend('GP3','zyf','RDG','zwl');
xlabel('zeta');
ylabel('the value of gap');
title('gap');
figure(1)
plot(zeta,performace_optimal_value,zeta,performance_GP3,zeta,performance_SP,zeta,performance_zwl,zeta,performance_zyf);
legend('optimal','GP3','RDG','zwl','zyf');
xlabel('zeta');
ylabel('the value of gap');
title('gap');
