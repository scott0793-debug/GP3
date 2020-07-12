clear ;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap

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
mu=[360,240,360,300,240,240,240,240,120,360,120,240,300,300,240,120,180,120,...
    120,180,60,300,300,60,180,180,300,360,240,480,360,300,360,240,240,360,...
    180,180,240,240,300,240,360,300,180,180,300,240,120,180,480,120,120,120,...
    180,240,180,120,240,240,240,360,300,360,120,180,180,300,120,240,240,240,...
    120,240,180,120]';
diag_sigma=[4,6,4,5,6,6,6,6,8,4,8,6,5,5,6,8,7,8,8,7,9,5,5,9,7,7,5,4,6,2,4,...
    5,4,6,6,4,7,7,6,6,5,6,4,5,7,7,5,6,8,7,2,8,8,8,7,6,7,8,6,6,6,4,5,4,8,7,...
    7,5,8,6,6,6,8,6,7,8];
covSigma=zeros(76);
for i=1:76
    covSigma(i,i)=diag_sigma(i);
end
b=zeros(24,1)';
b(1)=1;

number_zeta=100;
zeta=zeros(number_zeta,1);
for i=1:number_zeta
    zeta(i)=i;
end% zeta is tolerance
gap = 0.01;   
max_iter=400;
average_performance_zyf=zeros(23,number_zeta);
average_performance_zwl=zeros(23,number_zeta);
average_performance_SP=zeros(23,number_zeta);
average_performance_GP3=zeros(23,number_zeta);
average_iter_zyf=zeros(23,number_zeta);
average_iter_zwl=zeros(23,number_zeta);
average_iter_SP=zeros(23,number_zeta);
average_iter_GP3=zeros(23,number_zeta);
average_gap_zyf=zeros(23,number_zeta);
average_gap_zwl=zeros(23,number_zeta);
average_gap_RDG=zeros(23,number_zeta);
average_gap_GP3=zeros(23,number_zeta);
average_t_zyf=zeros(23,number_zeta);
average_t_zwl=zeros(23,number_zeta);
average_t_SP=zeros(23,number_zeta);
average_t_GP3=zeros(23,number_zeta);
average_number=10;
%% 终点从node2到node24
number_end_node=23;
for end_node=2:(number_end_node+1)
    b=zeros(24,1)';
    b(1)=1;
    b(end_node)=-1;
    %% zeta 从1到100
    for i=1:number_zeta
        i
        t_zyf_total=0;
        t_zwl_total=0;
        t_SP_total=0;
        t_GP3_total=0;
        iter_zyf_total=0;
        iter_zwl_total=0;
        iter_SP_total=0;
        iter_GP3_total=0;
        gap_zyf_total=0;
        gap_zwl_total=0;
        gap_RDG_total=0;
        gap_GP3_total=0;
        performance_zyf_total=0;
        performance_zwl_total=0;
        performance_SP_total=0;
        performance_GP3_total=0;
        %% 取平均值
        for z=1:average_number
            tic;
            [path_rsp_zyf,iter_zyf,gap_zyf,upper_bound_total_zyf,lower_bound_total_zyf,gap_total_zyf,gap_real_total_zyf,iteration_zyf]=func_rsp_zyf(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_zyf=toc;
            t_zyf_total=t_zyf+t_zyf_total;
            iter_zyf_total=iter_zyf_total+iter_zyf;
            gap_zyf_total=gap_zyf_total+gap_zyf;
            performance_zyf=mu'*path_rsp_zyf+zeta(i)*sqrt(path_rsp_zyf'*covSigma*path_rsp_zyf);
            performance_zyf_total=performance_zyf_total+performance_zyf;
            tic;
            [path_rsp_zwl,iter_zwl,gap_zwl,UB_zwl,LB_zwl,gap_total_zwf,gap_LB_UB_real,iteration_zwl]=func_rsp_zwl(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_zwl=toc;
            t_zwl_total=t_zwl_total+t_zwl;
            iter_zwl_total=iter_zwl+iter_zwl_total;
            gap_zwl_total=gap_zwl_total+gap_zwl;
            performance_zwl=mu'*path_rsp_zwl+zeta(i)*sqrt(path_rsp_zwl'*covSigma*path_rsp_zwl);
            performance_zwl_total=performance_zwl_total+performance_zwl;
            tic;
            [path_rsp_SP,iter_SP,gap_RDG,upper_bound_total_SP,lower_bound_total_SP,gap_RDG_total_SP,gap_RDG_real_SP,iteration_SP]=func_rsp_SP(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_SP=toc;
            t_SP_total=t_SP_total+t_SP;
            iter_SP_total=iter_SP_total+iter_SP;
            gap_RDG_total=gap_RDG_total+gap_RDG;
            performance_SP=mu'*path_rsp_SP+zeta(i)*sqrt(path_rsp_SP'*covSigma*path_rsp_SP);
            performance_SP_total=performance_SP_total+performance_SP;
            tic;
            [path_rsp_GP3,iter_GP3,gap_GP3,UB_GP3,LB_GP3,gap_total_GP3,gap_total_real_GP3,iteration_GP3]=func_GP3(A,b,mu,covSigma,zeta(i),gap);
            t_GP3=toc;
            t_GP3_total=t_GP3_total+t_GP3;
            iter_GP3_total=iter_GP3_total+iter_GP3;
            gap_GP3_total=gap_GP3_total+gap_GP3;
            performance_GP3=mu'*path_rsp_GP3+zeta(i)*sqrt(path_rsp_GP3'*covSigma*path_rsp_GP3);
            performance_GP3_total=performance_GP3_total+performance_GP3;
            table_GP3=table(iteration_GP3',UB_GP3',LB_GP3',gap_total_GP3',gap_total_real_GP3','VariableNames',{'iteration','UB','LB','gap','real_gap'});
            table_SP=table(iteration_SP',upper_bound_total_SP',lower_bound_total_SP',gap_RDG_total_SP',gap_RDG_real_SP','VariableNames',{'iteration','UB','LB','gap','real_gap'});
            table_zwl=table(iteration_zwl',UB_zwl',LB_zwl',gap_total_zwf',gap_LB_UB_real','VariableNames',{'iteration','UB','LB','gap','real_gap'});
            table_zyf=table(iteration_zyf',upper_bound_total_zyf',lower_bound_total_zyf',gap_total_zyf',gap_real_total_zyf','VariableNames',{'iteration','UB','LB','gap','real_gap'});
    
        end
        average_t_zyf(end_node-1,i)=t_zyf_total/average_number;
        average_t_zwl(end_node-1,i)=t_zwl_total/average_number;
        average_t_SP(end_node-1,i)=t_SP_total/average_number;
        average_t_GP3(end_node-1,i)=t_GP3_total/average_number;
        average_iter_zyf(end_node-1,i)=iter_zyf_total/average_number;
        average_iter_zwl(end_node-1,i)=t_zwl_total/average_number;
        average_iter_SP(end_node-1,i)=t_SP_total/average_number;
        average_iter_GP3(end_node-1,i)=iter_GP3_total/average_number;
        average_gap_zyf(end_node-1,i)=gap_zyf_total/average_number;
        average_gap_zwl(end_node-1,i)=gap_zwl_total/average_number;
        average_gap_RDG(end_node-1,i)=gap_RDG_total/average_number;
        average_gap_GP3(end_node-1,i)=gap_GP3_total/average_number;
        average_performance_zyf(end_node-1,i)=performance_zyf_total/average_number;
        average_performance_zwl(end_node-1,i)=performance_zwl_total/average_number;
        average_performance_SP(end_node-1,i)=performance_SP_total/average_number;
        average_performance_GP3(end_node-1,i)=performance_GP3_total/average_number;
    end
end
% figure(1)
% plot(zeta,gap_GP3,zeta,gap_zyf,zeta,gap_RDG,zeta,gap_zwl);
% legend('GP3','zyf','RDG','zwl');
% xlabel('zeta');
% ylabel('the value of gap');
% title('gap');
% figure(2)
% plot(zeta,performance_GP3,zeta,performance_zyf,zeta,performance_SP,zeta,performance_zwl);
% legend('GP3','zyf','RDG','zwl');
% xlabel('zeta');
% ylabel('the value of performance');
% title('performance');
% figure(3)
% plot(zeta,t_GP3,zeta,t_zyf,zeta,t_SP,zeta,t_zwl);
% legend('GP3','zyf','SP','zwl');
% xlabel('zeta');
% ylabel('time');
% title('time');
for i=1:number_end_node
    figure(i);
    plot(zeta,average_performance_GP3(i,:),zeta,average_performance_SP(i,:),zeta,average_performance_zwl(i,:),zeta,average_performance_zyf(i,:));
    legend('GP3','SP','zwl','zyf');
    xlabel('zeta');
    ylabel('the value of performance')
end
for i=(number_end_node+1):(2*number_end_node)
    figure(i);
    plot(zeta,average_t_GP3((i-number_end_node),:),zeta,average_t_SP((i-number_end_node),:),zeta,average_t_zwl((i-number_end_node),:),zeta,average_t_zyf((i-number_end_node),:));
    legend('GP3','SP','zwl','zyf');
    xlabel('zeta');
    ylabel('time')
end
for i=(2*number_end_node+1):(3*number_end_node)
    figure(i);
    plot(zeta,average_iter_GP3((i-2*number_end_node),:),zeta,average_iter_SP((i-2*number_end_node),:),zeta,average_iter_zwl((i-2*number_end_node),:),zeta,average_iter_zyf((i-2*number_end_node),:));
    legend('GP3','SP','zwl','zyf');
    xlabel('zeta');
    ylabel('itertation')
end
for i=(3*number_end_node+1):(4*number_end_node)
    figure(i);
    plot(zeta,average_gap_GP3((i-3*number_end_node),:),zeta,average_gap_RDG((i-3*number_end_node),:),zeta,average_gap_zwl((i-3*number_end_node),:),zeta,average_gap_zyf((i-3*number_end_node),:));
    legend('GP3','SP','zwl','zyf');
    xlabel('zeta');
    ylabel('itertation')
end