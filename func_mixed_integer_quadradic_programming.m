function [path_rsp_MIQP]=func_mixed_integer_quadradic_programming(A,b,covSigma,mu,gap,zeta)

%% calculate #nodes and #edges out of the map
[num_nodes,num_edges]=size(A); %num_nodes is the number of nodes.num_edges is the number of edges

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


%% initialize t=(t_min+t_max)/2, and calculate matrix H, lambda_min and convert to Q and l;
%% initialize x_rsp as x_let

H = zeta^2*covSigma-mu*mu';

lambda=eig(H);
lambda_min = min(lambda);
Q = H-(lambda_min-1)*eye(num_edges);
%l = 2*t*mu+(lambda_min-1)*ones(num_edges,1); %no need to initialize l
%here, because we will re-initialize l in the while loop

path_rsp_MIQP = x_let;
while (t_max-t_min)>gap 
    t=(t_min+t_max)/2;
    % update l, 
    l = 2*t*mu+(lambda_min-1)*ones(num_edges,1);  
    % find the optimal path for the constrained quadratic program
    A_mu=zeros(num_nodes+1,num_edges);
    A_mu(1:num_nodes,:)=A;
    A_mu(num_nodes+1,:)=mu';
    lower_bound_b=zeros(num_nodes+1,1);
    upper_bound_b=zeros(num_nodes+1,1);
    lower_bound_b(1:num_nodes,1)=b;
    upper_bound_b(1:num_nodes,1)=b;
    lower_bound_b(num_nodes,1)=0;
    upper_bound_b(num_nodes,1)=t;
    
    [~,x_rsp_propose] = func_rsp_path(Q,l,A,b,mu,t);
    rsp_min=mu'*x_rsp_propose+zeta*sqrt(x_rsp_propose'*covSigma*x_rsp_propose);
    if rsp_min<t
        t_max = rsp_min;
        path_rsp_MIQP = x_rsp_propose;
    else
        t_min = t;

    gap_GP3=(t_max-t_min);

    iter_GP3=iter_GP3+1;
end
end