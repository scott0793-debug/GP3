function [obj_min,x_rsp]=func_simplepathleastvar(A,b,mu,covSigma)

%% prepare the parameters, and solve the continous relaxation problem
[~,num_edges]=size(A); %num_nodes is the number of nodes.num_edges is the number of edges
lb=zeros(num_edges,1);
ub=ones(num_edges,1);
x_rsp_quadprog=quadprog(2*covSigma,[],[],[],A,b,lb,ub); % x_rsp_quadprog is the optimal solution of the (continously) relaxed problem




%% find the set of elementary paths: path_
mu_dynamic = mu;
counter=0;
while ones(num_edges,1)'*x_rsp_quadprog>0.1 % if there is remaining x_rsp_quadprog left, continue;
    counter=counter+1;
    x_rsp_quadprog_value=x_rsp_quadprog;
    for i=1:num_edges
        
        if x_rsp_quadprog(i)<=0.001
            x_rsp_quadprog_value(i)=100;
            mu_dynamic(i)=100000;
        end         
    end
    clear i
    [edge_value,edge_index] = min(x_rsp_quadprog_value);  % find the edge list with the smallest value of x_rsp_quadprog  


    path_dijkstraPP=func_dijkstraPP(A,b,mu_dynamic,edge_index);% 找到最短路径
    
    path_p_i=zeros(num_edges,1);%initialize path_p_i as the elementary path
    length_path_dijkstraPP = length(path_dijkstraPP);
    for i=1:length_path_dijkstraPP
        path_p_i(path_dijkstraPP(i))=1;% set path_p_i according to the previously calculated path from func_DijkstraPP  %% too complex
    end 
    clear i
    path_element_set(:,counter)=path_p_i;% save the elementary paths
    
    x_rsp_quadprog=x_rsp_quadprog-edge_value*path_p_i;% update x_rsp_quadprog; decrease the original x_rsp_quadprog value
    
end


 


%% select the best elementary path as the solution
objective_values=zeros(counter,1);
for i=1:counter 
    objective_values(i)=path_element_set(:,i)'*covSigma*path_element_set(:,i);% calculate the objective value which is x'*Q*x+l'*x
end
[obj_min,index] = min(objective_values);
%position=find(objective_values==obj_min);
x_rsp=path_element_set(:,index);% 找到方差最小的路径





end
