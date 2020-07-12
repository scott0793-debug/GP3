function [path_rsp_CG,iter_CG,gap_RIG,gap_RDG]=func_rsp_CG(A,b,mu,covSigma,zeta,max_iter,gap)
%% extract information a map, and perform cholesky decomposition over Sigma
[~,num_edges]=size(A); %num_edges is the number of edges

P_low_triangular=chol(covSigma)';% cholesky decomposition; get the lower triangular matrix

%% initialization ; f 是path时间的均值和方差的线性组合值 f_D是(mu+P*lambda)'*X的最优值
lambda=zeros(num_edges,1);
for i=1:num_edges
    lambda(i)=0.001;
end
upper_fD=100000;% upper bound for fD
upper_f=10000;% upper_f is  upper bound for f*
lower_fD=-10000;% low bound for fD
iter_CG=1;
gap_RIG=1;%relative duality gap
gap_RDG=1;% relative inner gap
x_let_dijkstra=func_dijkstra(A,b,mu);%a standard shorest path
path_rsp_CG=zeros(num_edges,1);
length_x_let_dijkstra=length(x_let_dijkstra);
for i=1:length_x_let_dijkstra
    path_rsp_CG(x_let_dijkstra(i))=1;
end
%% constraint generation
while iter_CG<max_iter && gap_RIG>=gap && gap_RDG>=gap
    dynamic_mu=mu+P_low_triangular*lambda;
    x_let_dijkstra=func_dijkstra(A,b,dynamic_mu);%a standard shorest path
    dynamic_path_rsp_CG=zeros(num_edges,1);
    length_x_let_dijkstra=length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        dynamic_path_rsp_CG(x_let_dijkstra(:,i))=1;
    end
    % update low bound of fD
    optimal_fD=dynamic_mu'*dynamic_path_rsp_CG;
    lower_fD=max(lower_fD,optimal_fD);
    % update upper bound of f
    optimal_f=mu'*dynamic_path_rsp_CG+zeta*sqrt(dynamic_path_rsp_CG'*covSigma*dynamic_path_rsp_CG);
    if upper_f>optimal_f
        upper_f=optimal_f;
        path_rsp_CG=dynamic_path_rsp_CG;
    end
 %% solve the master problem
    cvx_begin

            variables lambda(num_edges) ;
            variables value_gamma

            maximize( value_gamma );

            subject to
                value_gamma<=mu'*dynamic_path_rsp_CG+dynamic_path_rsp_CG'*P_low_triangular*lambda;
                lambda'*lambda<=zeta^2;


    cvx_end
    % update upper bound of fD

    upper_fD=min(upper_fD,value_gamma);
    gap_RDG=(upper_f-lower_fD)/lower_fD;
    gap_RIG=(upper_fD-lower_fD)/lower_fD;
    iter_CG=iter_CG+1;
end
end



