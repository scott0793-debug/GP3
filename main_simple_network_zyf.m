clear ;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap
A=[1,-1,0,0,0,0;1,0,-1,0,0,0;0,1,0,-1,0,0;0,1,0,0,-1,0;0,-1,1,0,0,0;0,0,1,-1,0,0;0,0,0,1,0,-1;0,0,0,-1,1,0;0,0,0,0,1,-1]';
b=[1,0,0,0,0,-1]';
mu=[5,4,1,2,2,2,6,4,4]';
covSigma=[5,1.2,0.75,0.75,-0.325,0.75,0.5,0.7625,1.5;1.2,3,0.875,-0.5,0.55,-0.8,-1.125,0.95,0.7875;...
    0.75,0.875,2,1,1,0.8,-0.3,0.3875,0.5125;0.75,-0.5,1,3,1.25,0.9,0.5,-0.325,1;...
    -0.325,0.55,1,1.25,2,0.7,-0.5,0.375,0.275;0.75,-0.8,0.8,0.9,0.7,3.5,0.625,-0.3625,-0.1375;...
    0.5,-1.125,-0.3,0.5,-0.5,0.625,2,0.8125,1;0.7625,0.95,0.3875,-0.325,0.375,-0.3625,0.8125,3.5,-0.6125;...
    1.5,0.7875,0.5125,1,0.275,-0.1375,1,-0.6125,4];
zeta = 1;      % zeta is tolerance
gap = 0.01;   
max_iter=100;
tic;
[path_rsp_zyf,iter_zyf,gap_zyf,upper_bound_total_zyf,lower_bound_total_zyf,gap_total_zyf,gap_real_total_zyf,iteration_zyf]=func_rsp_zyf(A,b,mu,covSigma,zeta,max_iter,gap);
t_zyf=toc;
performance_zyf=mu'*path_rsp_zyf+zeta*sqrt(path_rsp_zyf'*covSigma*path_rsp_zyf);
tic;
[path_rsp_zwl,iter_zwl,gap_zwl,UB_zwl,LB_zwl,gap_total_zwf,gap_LB_UB_real,iteration_zwl]=func_rsp_zwl_Lw(A,b,mu,covSigma,zeta,max_iter,gap);
t_zwl=toc;
performance_zwl=mu'*path_rsp_zwl+zeta*sqrt(path_rsp_zwl'*covSigma*path_rsp_zwl);
tic;
[path_rsp_SP,iter_SP,gap_RDG,upper_bound_total_SP,lower_bound_total_SP,gap_RDG_total_SP,gap_RDG_real_SP,iteration_SP]=func_rsp_SP(A,b,mu,covSigma,zeta,max_iter,gap);
t_SP=toc;
performance_SP=mu'*path_rsp_SP+zeta*sqrt(path_rsp_SP'*covSigma*path_rsp_SP);
tic;
[path_rsp_GP3,iter_GP3,gap_GP3,UB_GP3,LB_GP3,gap_total_GP3,gap_total_real_GP3,iteration_GP3]=func_GP3(A,b,mu,covSigma,zeta,gap);
t_GP3=toc;
performance_GP3=mu'*path_rsp_GP3+zeta*sqrt(path_rsp_GP3'*covSigma*path_rsp_GP3);
table_GP3=table(iteration_GP3',UB_GP3',LB_GP3',gap_total_GP3',gap_total_real_GP3','VariableNames',{'iteration','UB','LB','gap','real_gap'});
table_SP=table(iteration_SP',upper_bound_total_SP',lower_bound_total_SP',gap_RDG_total_SP',gap_RDG_real_SP','VariableNames',{'iteration','UB','LB','gap','real_gap'});
table_zwl=table(iteration_zwl',UB_zwl',LB_zwl',gap_total_zwf',gap_LB_UB_real','VariableNames',{'iteration','UB','LB','gap','real_gap'});
table_zyf=table(iteration_zyf',upper_bound_total_zyf',lower_bound_total_zyf',gap_total_zyf',gap_real_total_zyf','VariableNames',{'iteration','UB','LB','gap','real_gap'});
