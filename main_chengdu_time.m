clear ;
clc;
close all;
%% 输入的内容，包括1.地图A,b 2.交通特性 mu ,Sigma 3.用户特性 zeta 4.算法终止条件num，gap
load('holiday_off_peak_covSigma.mat')
load('holiday_off_peak_mu.mat')
load('chengdu_mapA.mat');
zeta=1;
gap=0.1;
max_iter=5;
counter=0;
[num_nodes,num_edges]=size(A);
lambda01=1;
lambda02=ones(1,num_edges);
lambda03=ones(1,num_edges);
omega=ones(1,num_nodes);
H = zeta^2*covSigma-mu*mu';
t=1;
lambda=eig(H);
lambda_min = min(lambda);
Q = H-(lambda_min-1)*eye(num_edges);
inv_Q=-inv(Q);
t1=0;
t2=0;
t3=0;
t4=0;

for iter=1:11
    b=zeros(num_nodes,1);
    total_nodes=[1:1:num_nodes];
    OD=total_nodes(randperm(length(total_nodes),2));
    b(min(OD))=1;
    b(max(OD))=-1;

    % step 1: calculate the LET path: x_let;
    x_let_dijkstra = func_dijkstra(A,b,mu);% A: map, mu: travel time; b: OD
    % convert to LET path representation form
    x_let=zeros(num_edges,1); %x_let is the vector representation of the LET path
    length_x_let_dijkstra = length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        x_let(x_let_dijkstra(i))=1;
    end

    % step 2: calculate the smallest variance path: x_var, and the variance
    % value (var_min)
    [~,x_var]=func_simplepathleastvar(A,b,mu,covSigma);

    % step 3: calculate t_min and t_max
    t_min=mu'*x_let+sqrt(x_var'*covSigma*x_var); % t_min is the sum of the two min values
    t_max=mu'*x_let+sqrt(x_let'*covSigma*x_let); % t_max is one of the feasible path's value which the LET path
    t=(t_min+t_max)/2;
    l = 2*t*mu+(lambda_min-1)*ones(num_edges,1); 
    tic;
    x1=func_GP3_accelerate(inv_Q,lambda_min,lambda01,mu,lambda02,lambda03,omega,A,t_min,t_max,b,num_edges,num_nodes,covSigma);
    t1=toc+t1;
end
res1=t1/11;