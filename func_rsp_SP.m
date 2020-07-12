function [path_rsp_SP,iter_SP,gap_RDG,upper_bound_total_SP,lower_bound_total_SP,gap_RDG_total_SP,gap_RDG_real_SP,iteration_SP]=func_rsp_SP(A,b,mu,covSigma,zeta,max_iter,gap)
%% extract information a map, and perform cholesky decomposition over Sigma
[~,num_edges]=size(A); %num_edges is the number of edges

P_low_triangular=chol(covSigma)';% cholesky decomposition; get the lower triangular matrix

%% initialization ; f 是path时间的均值和方差的线性组合值 f_D是(mu+P*lambda)'*X的最优值
lambda=zeros(num_edges,1);
for i=1:num_edges
    lambda(i)=zeta/sqrt(num_edges);
end

upper_f=1000000;% upper_f is  upper bound for f*
lower_fD=-1000000;% low bound for fD
iter_SP=1;

gap_RDG=1;% relative inner gap
x_let_dijkstra=func_dijkstra(A,b,mu);%a standard shorest path
path_rsp_SP=zeros(num_edges,1);
length_x_let_dijkstra=length(x_let_dijkstra);
for i=1:length_x_let_dijkstra
    path_rsp_SP(x_let_dijkstra(i))=1;
end
%% constraint generation
while iter_SP<max_iter &&  (upper_f-lower_fD)>=gap
    dynamic_mu=mu+P_low_triangular*lambda;
    x_let_dijkstra=func_dijkstra(A,b,dynamic_mu);%a standard shorest path
    dynamic_path_rsp_SP(:,iter_SP)=zeros(num_edges,1);
    length_x_let_dijkstra=length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        dynamic_path_rsp_SP(x_let_dijkstra(:,i),iter_SP)=1;
    end
    %  compute y and update lower bound of fD
    optimal_fD=dynamic_mu'*dynamic_path_rsp_SP(:,iter_SP);
   
    y=lambda'*P_low_triangular'*dynamic_path_rsp_SP(:,iter_SP)*lambda/(zeta^2);
    lower_fD=max(lower_fD,optimal_fD);
    % update upper bound of f
    optimal_f=mu'*dynamic_path_rsp_SP(:,iter_SP)+zeta*sqrt(dynamic_path_rsp_SP(:,iter_SP)'*covSigma*dynamic_path_rsp_SP(:,iter_SP));
    if upper_f>optimal_f
        upper_f=optimal_f;
        path_rsp_SP=dynamic_path_rsp_SP(:,iter_SP);
    end
 %% update lambda
    for j=1:iter_SP
        dual_value(j)=mu'*dynamic_path_rsp_SP(:,j)+zeta*sqrt(dynamic_path_rsp_SP(:,j)'*covSigma*dynamic_path_rsp_SP(:,j));
    end
    approximation_optimal_dual_value=min(dual_value);
    delta=0.1;
    eta_numerator=delta*(approximation_optimal_dual_value-optimal_fD);
    eta_denominator=(y-P_low_triangular'*dynamic_path_rsp_SP(:,iter_SP))'*(y-P_low_triangular'*dynamic_path_rsp_SP(:,iter_SP));
    eta=eta_numerator/eta_denominator;
    sub_lambda=lambda+eta*(P_low_triangular'*dynamic_path_rsp_SP(:,iter_SP)-y);
    value_sub_lambda=sqrt(sub_lambda'*sub_lambda);
    lambda=zeta*sub_lambda/value_sub_lambda;
    % update gap_RDG

    gap_RDG=(upper_f-lower_fD)/upper_f;
    gap_RDG_total_SP(iter_SP)=upper_f-lower_fD;
    gap_RDG_real_SP(iter_SP)=upper_f-lower_fD;
    upper_bound_total_SP(iter_SP)=upper_f;
    lower_bound_total_SP(iter_SP)=lower_fD;
    iteration_SP(iter_SP)=iter_SP;
    
    iter_SP=iter_SP+1;
end
end



