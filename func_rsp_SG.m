function [path_rsp_SG,iter_SG,gap_SG]=func_rsp_SG(A,b,mu,Sigma,zeta,num,gap)
%relative duality gap threshold epsilon,num is maximum number of iteration 
%% step 1,variance-covariance maxtrix decomposition
[~,num_edges]=size(A); %num_edges is the number of edges
[eig_Sigma,eig_Lambda] = eig(Sigma);
%% step2 initialization
iter_SG=0;%set iteraton number
lower_bound_Z=-1000; %initial lower bounder of Z
Lagrangian_multiplier(:,iter_SG+1)=ones(1,num_edges);%initial lagrangian multiplier
x_let_dijkstra=func_dijkstra(A,b,mu);%a standard shorest path
x_let=zeros(num_edges,1);
length_x_let_dijkstra=length(x_let_dijkstra);
for i=1:length_x_let_dijkstra
    x_let(x_let_dijkstra(:,i))=1;
end
upper_bound_Z=mu'*x_let+zeta*sqrt(x_let'*Sigma*x_let);%compute Z as the upper bound
gap_SG=(upper_bound_Z-lower_bound_Z)/upper_bound_Z;
while iter_SG<num && gap_SG>gap
    %% step 3 solve lagtangian telaxation problems

    dynamic_mu=mu+eig_Sigma*Lagrangian_multiplier(:,iter_SG+1);
    dynamic_mu=abs(dynamic_mu);
    [x_let_dijkstra]=func_dijkstra(A,b,dynamic_mu);%a standard shorest path
    x_let=zeros(num_edges,1);
    length_x_let_dijkstra=length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        x_let(x_let_dijkstra(:,i))=1;
    end
    optimal_upper_bound_Z=mu'*x_let+zeta*sqrt(x_let'*Sigma*x_let);% update upper bound of Z
    Lagrangian_function_optimal_value=(mu+eig_Sigma*Lagrangian_multiplier(:,iter_SG+1))'*x_let;
    %% update upper_bound_Z and save x_let
    if upper_bound_Z>=optimal_upper_bound_Z
        upper_bound_Z=optimal_upper_bound_Z;
        x_SG_total(:,iter_SG+1)=x_let';
    else 
        if iter_SG==0
            x_SG_total(:,iter_SG+1)=x_let';
        else
           x_SG_total(:,iter_SG+1)=x_SG_total(:,iter_SG);
        end
    end
    Lagrangian_function_value=Lagrangian_function_optimal_value;
    %% update low_bound_Z
    if Lagrangian_function_value>=lower_bound_Z
        lower_bound_Z=Lagrangian_function_value;
    end
    % update lagrangian multiplier
    partial_coefficient_lagrangian_multiplier=eig_Sigma'*x_let;
    coefficient_lagrangian_multiplier=abs(upper_bound_Z-lower_bound_Z)/(partial_coefficient_lagrangian_multiplier'*partial_coefficient_lagrangian_multiplier);
    beta=0.1;
    Lagrangian_multiplier(:,iter_SG+2)=Lagrangian_multiplier(:,iter_SG+1)+beta*coefficient_lagrangian_multiplier*partial_coefficient_lagrangian_multiplier;
    for i=1:num_edges
        if Lagrangian_multiplier(i,iter_SG+2)>zeta*sqrt(eig_Lambda(i,i))
            Lagrangian_multiplier(i,iter_SG+2)=Lagrangian_multiplier(i,iter_SG+1);
        elseif Lagrangian_multiplier(i,iter_SG+2)<-zeta*sqrt(eig_Lambda(i,i))
            Lagrangian_multiplier(i,iter_SG+2)=Lagrangian_multiplier(i,iter_SG+1);
        end
    
    end
    iter_SG=iter_SG+1;
    path_rsp_SG=x_SG_total(:,end);
    gap_SG=(upper_bound_Z-lower_bound_Z)/upper_bound_Z;
end
end


