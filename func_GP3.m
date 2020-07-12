function [path_rsp_GP3,iter_GP3,gap_GP3,UB_GP3,LB_GP3,gap_total_GP3,gap_total_real_GP3,iteration_GP3]=func_GP3(A,b,mu,covSigma,zeta,gap)

%% calculate #nodes and #edges out of the map
[~,num_edges]=size(A); %num_nodes is the number of nodes.num_edges is the number of edges

%% calculate t_min and t_max (lower bound and upper bound)
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
t_min=mu'*x_let+zeta*sqrt(x_var'*covSigma*x_var); % t_min is the sum of the two min values
t_max=mu'*x_let+zeta*sqrt(x_let'*covSigma*x_let); % t_max is one of the feasible path's value which the LET path

gap_GP3=t_max-t_min;
UB_GP3=t_max;
LB_GP3=t_max;
gap_total_GP3=(t_max-t_min)/t_max;
gap_total_real_GP3=t_max-t_min;
%% initialize t=(t_min+t_max)/2, and calculate matrix H, lambda_min and convert to Q and l;
%% initialize x_rsp as x_let

H = zeta^2*covSigma-mu*mu';

lambda=eig(H);
lambda_min = min(lambda);
Q = H-(lambda_min-1)*eye(num_edges);
%l = 2*t*mu+(lambda_min-1)*ones(num_edges,1); %no need to initialize l
%here, because we will re-initialize l in the while loop

path_rsp_GP3 = x_let;
iter_GP3=1;
iteration_GP3=1;
%% loop and update t until convergence. if no feasible path, set t_min = t; if feasible path, set t_max = rsp_min;
while (t_max-t_min)>gap 
    t=(t_min+t_max)/2;
    % update l, 
    l = 2*t*mu+(lambda_min-1)*ones(num_edges,1);  
    % find the optimal path for the constrained quadratic program
    [~,x_rsp_propose] = func_rsp_path(Q,l,A,b,mu,t);
    %[~,x_rsp_propose] = func_rsp_path(Q,l,A,b,mu,t-sqrt(x_var'*covSigma*x_var));
    rsp_min=mu'*x_rsp_propose+zeta*sqrt(x_rsp_propose'*covSigma*x_rsp_propose);
    
    %% 20200617 I added the following logic so that whenever we get a better rsp, we store it into the GP3 solver
    if rsp_min < t_max
        t_max = rsp_min;
        path_rsp_GP3 = x_rsp_propose;
    end
    
    
    if rsp_min < t_min
        t_min = rsp_min;
    %% end of the previous edits. 
    elseif rsp_min<t
        t_max = rsp_min;
        path_rsp_GP3 = x_rsp_propose;
    else
        t_min = t;
    end
    iteration_GP3(iter_GP3)=iter_GP3;
    
    gap_GP3=(t_max-t_min);
    UB_GP3(iter_GP3)=t_max;
    LB_GP3(iter_GP3)=t_min;
    gap_total_GP3(iter_GP3)=(t_max-t_min)/t_max;
    gap_total_real_GP3(iter_GP3)=t_max-t_min;
    iteration_GP3(iter_GP3)=iter_GP3;
    iter_GP3=iter_GP3+1;
end


               

















