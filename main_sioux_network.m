clear ;
clc;
close all;

%% The input contents include 1. Map A, b 2. Traffic characteristics mu, sigma
%%3. User characteristics zeta 4. Algorithm termination conditions num, gap
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

load('mu.mat');
load('covSigma.mat');

number_zeta=60;
zeta=zeros(number_zeta,1);
for i=1:number_zeta
    zeta(i)=i;
end% zeta is tolerance
gap = 0.01;   
max_iter=200;
% average  value of linear combination of the path's expected travel time and standard deviation
average_performance_SG=zeros(23,number_zeta);
average_performance_LR=zeros(23,number_zeta);
average_performance_SP=zeros(23,number_zeta);
average_performance_GP3=zeros(23,number_zeta);
% average value of iteration
average_iter_SG=zeros(23,number_zeta);
average_iter_LR=zeros(23,number_zeta);
average_iter_SP=zeros(23,number_zeta);
average_iter_GP3=zeros(23,number_zeta);
% average value of relative gap
average_gap_SG=zeros(23,number_zeta);
average_gap_LR=zeros(23,number_zeta);
average_gap_RDG=zeros(23,number_zeta);
average_gap_GP3=zeros(23,number_zeta);
% average value of running time
average_t_SG=zeros(23,number_zeta);
average_t_LR=zeros(23,number_zeta);
average_t_SP=zeros(23,number_zeta);
average_t_GP3=zeros(23,number_zeta);

%% The end point is from node2 to node24
number_end_node=23;
for end_node=2:(1+number_end_node)
    b=zeros(24,1)';
    b(1)=1;
    b(end_node)=-1;
    path_rsp_GP3_total=zeros(76,100);
    %%zeta goes from 1 to 100
    for i=1:number_zeta
        %initial the value of four algorithms' total running time
        t_SG_total=0;
        t_LR_total=0;
        t_SP_total=0;
        t_GP3_total=0;
        %initial the value of four algorithms' total iterations
        iter_SG_total=0;
        iter_LR_total=0;
        iter_SP_total=0;
        iter_GP3_total=0;
        %initial the value of four algorithms' total gap
        gap_SG_total=0;
        gap_LR_total=0;
        gap_RDG_total=0;
        gap_GP3_total=0;
        %initial the value of four algorithms' mean-std RSP
        performance_SG_total=0;
        performance_LR_total=0;
        performance_SP_total=0;
        performance_GP3_total=0;
        %% running average_number times and then  get the average value of all above value
        average_number=1000;
        for z=1:average_number
            tic;
            [path_rsp_SG,iter_SG,gap_SG]=func_rsp_SG(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_SG=toc;
            t_SG_total=t_SG+t_SG_total;
            iter_SG_total=iter_SG_total+iter_SG;
            gap_SG_total=gap_SG_total+gap_SG;
            performance_SG=mu'*path_rsp_SG+zeta(i)*sqrt(path_rsp_SG'*covSigma*path_rsp_SG);
            performance_SG_total=performance_SG_total+performance_SG;
            tic;
            [path_rsp_LR,iter_LR,gap_LR] = func_rsp_LR(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_LR=toc;
            t_LR_total=t_LR_total+t_LR;
            iter_LR_total=iter_LR+iter_LR_total;
            gap_LR_total=gap_LR_total+gap_LR;
            performance_LR=mu'*path_rsp_LR+zeta(i)*sqrt(path_rsp_LR'*covSigma*path_rsp_LR);
            performance_LR_total=performance_LR_total+performance_LR;
            tic;
            [path_rsp_SP,iter_SP,gap_RDG]=func_rsp_SP(A,b,mu,covSigma,zeta(i),max_iter,gap);
            t_SP=toc;
            t_SP_total=t_SP_total+t_SP;
            iter_SP_total=iter_SP_total+iter_SP;
            gap_RDG_total=gap_RDG_total+gap_RDG;
            performance_SP=mu'*path_rsp_SP+zeta(i)*sqrt(path_rsp_SP'*covSigma*path_rsp_SP);
            performance_SP_total=performance_SP_total+performance_SP;
            tic;
            [path_rsp_GP3,iter_GP3,gap_GP3]=func_GP3(A,b,mu,covSigma,zeta(i),gap);
            path_rsp_GP3_total(:,i)=path_rsp_GP3;
            t_GP3=toc;
            t_GP3_total=t_GP3_total+t_GP3;
            iter_GP3_total=iter_GP3_total+iter_GP3;
            gap_GP3_total=gap_GP3_total+gap_GP3;
            performance_GP3=mu'*path_rsp_GP3+zeta(i)*sqrt(path_rsp_GP3'*covSigma*path_rsp_GP3);
            performance_GP3_total=performance_GP3_total+performance_GP3;

        end
        average_t_SG(end_node-1,i)=t_SG_total/average_number;
        average_t_LR(end_node-1,i)=t_LR_total/average_number;
        average_t_SP(end_node-1,i)=t_SP_total/average_number;
        average_t_GP3(end_node-1,i)=t_GP3_total/average_number;
        average_iter_SG(end_node-1,i)=iter_SG_total/average_number;
        average_iter_LR(end_node-1,i)=iter_LR_total/average_number;
        average_iter_SP(end_node-1,i)=iter_SP_total/average_number;
        average_iter_GP3(end_node-1,i)=iter_GP3_total/average_number;
        average_gap_SG(end_node-1,i)=gap_SG_total/average_number;
        average_gap_LR(end_node-1,i)=gap_LR_total/average_number;
        average_gap_RDG(end_node-1,i)=gap_RDG_total/average_number;
        average_gap_GP3(end_node-1,i)=gap_GP3_total/average_number;
        average_performance_SG(end_node-1,i)=performance_SG_total/average_number;
        average_performance_LR(end_node-1,i)=performance_LR_total/average_number;
        average_performance_SP(end_node-1,i)=performance_SP_total/average_number;
        average_performance_GP3(end_node-1,i)=performance_GP3_total/average_number;
    end
end
