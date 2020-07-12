function [obj_min,x_rsp]=func_rsp_path(Q,l,A,b,mu,t)

small_value = 0.00005;
MAX_objective_value=10000000;
MAX_edge_value=10000;
%% prepare the parameters, and solve the continous relaxation problem
[~,num_edges]=size(A); %num_nodes is the number of nodes.num_edges is the number of edges
lb=zeros(num_edges,1)+small_value*ones(num_edges,1);
ub=ones(num_edges,1)-small_value*ones(num_edges,1);
Q=(Q+Q')/2;
x_rsp_quadprog=quadprog(2*Q,l,mu',t,A,b,lb,ub); % x_rsp_quadprog is the optimal solution of the (continously) relaxed problem
%x_rsp_quadprog=quadprog(2*Q,l,[],[],A,b,lb,ub); % x_rsp_quadprog is the optimal solution of the (continously) relaxed problem


dynamic_mu=mu;
%% find the set of elementary paths: path_
counter=0;
while ones(num_edges,1)'*x_rsp_quadprog>0.01 % if there is remaining x_rsp_quadprog left, continue;
    counter=counter+1;
    
    length_x_rsp_quadprog=length(x_rsp_quadprog);
    dynamic_x_rsp_quadprog=x_rsp_quadprog;
    for i=1:length_x_rsp_quadprog
        if x_rsp_quadprog(i)<0.0001
            dynamic_mu(i)=MAX_edge_value;
            dynamic_x_rsp_quadprog(i)=100;
        end
    end
    [edge_value,edge_index] = min(dynamic_x_rsp_quadprog);  % find the edge list with the smallest value of x_rsp_quadprog  
    %path_dijkstraPP=func_dijkstraPP(A,b,dynamic_mu,edge_index);% 找到最短路径
    path_dijkstraPP=func_dijkstraPP(A,b,dynamic_mu,edge_index);% 找到mu+zeta*sqrt(sigma)路径
    path_p_i=zeros(num_edges,1);%initialize path_p_i as the elementary path
    lenght_path_dijkstraPP=length(path_dijkstraPP);
    for i=1:lenght_path_dijkstraPP
        path_p_i(path_dijkstraPP(i))=1;% set path_p_i according to the previously calculated path from func_DijkstraPP
    end
    
    x_rsp_quadprog=x_rsp_quadprog-edge_value*path_p_i;% update x_rsp_quadprog
    lenght_x_rsp_quadprog=length(x_rsp_quadprog);
    for i=1:lenght_x_rsp_quadprog
        if x_rsp_quadprog(i)<0.001
            dynamic_mu(i)=MAX_edge_value;
        end
    end
    path_element_set(:,counter)=path_p_i;% save the elementary paths
end


 


%% select the best elementary path as the solution
objective_values=zeros(counter,1);
for i=1:counter 
    path_mean_value=path_element_set(:,i)'*mu;
    if path_mean_value>t
        objective_values(i)=MAX_objective_value;
    else
        objective_values(i)=path_element_set(:,i)'*Q*path_element_set(:,i)+l'*path_element_set(:,i); % calculate the objective value which is x'*Q*x+l'*x
    end
end
[obj_min,index] = min(objective_values);
%position=find(objective_values==obj_min);
x_rsp=path_element_set(:,index);% 找到方差最小的路径





end
