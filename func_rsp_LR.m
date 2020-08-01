function [path_rsp_LR,iter_LR,gap_LR] = func_rsp_LR(A,b,mu,covSigma,zeta,max_iter,gap)
%% extract information a map, and perform cholesky decomposition over Sigma
[~,num_edges]=size(A); %num_edges is the number of edges

L_low_triangular=chol(covSigma)';% cholesky decomposition; get the lower triangular matrix


%% initialization, calculate lowerbound and upperbound of the problem
lambda=rand(num_edges,1);%randomly generate the initial values for lamda,v (0,1)
nu=rand;
iter_LR=1;% initialize iteration number
x_let_dijkstra=func_dijkstra(A,b,mu);%a standard shorest path
x_let=zeros(num_edges,1);
length_x_let_dijkstra=length(x_let_dijkstra);
for i=1:length_x_let_dijkstra
    x_let(x_let_dijkstra(:,i))=1;
end
y_pri=(L_low_triangular'*x_let)'*(L_low_triangular'*x_let);%the travel time variance of the shortest path
 % y' is the travel time variance of the shortest distance path
upperbound=mu'*x_let+zeta*sqrt(x_let'*covSigma*x_let);%calculate UB using Eq.40
lowerbound=-100;
gap_LR=(upperbound-lowerbound)/upperbound;
L_high_triangular=L_low_triangular';
%% step2 slove the decomposed dual problems
while iter_LR<max_iter && abs(gap_LR)>gap
    %% solve a transformed shortest path problem

    lambda_L_high_triangular=lambda'*L_high_triangular;
    length_lambda_L_high_triangular=length(lambda_L_high_triangular);
    nan_position=isnan(lambda_L_high_triangular);
    length_nan_position=length(nan_position);
    for i=1:length_nan_position
        if nan_position(i)~=0
            lambda_L_high_triangular(i)=0;
            lambda(i)=0;
        end
    end
    for i=1:length_lambda_L_high_triangular
        if lambda_L_high_triangular(i)<0
            lambda_L_high_triangular(i)=-lambda_L_high_triangular(i);
            
        end
    end

    dynamic_mu=mu+lambda_L_high_triangular';

   
    x_let_dijkstra=func_dijkstra(A,b,dynamic_mu);%a standard shorest path using sub-function using a standard shortest-path algorithm
    x_let_dynamic=zeros(num_edges,1);
    length_x_let_dijkstra=length(x_let_dijkstra);
    for i=1:length_x_let_dijkstra
        x_let_dynamic(x_let_dijkstra(i))=1;
    end
    path_rsp_LR=x_let_dynamic;
    
    %% calculate Lx,Ly,Lw
    
    Lx=dynamic_mu'*x_let_dynamic;
    Lx_total(iter_LR)=Lx;
  
    omega=lambda/(2*nu);%  calculate omega ,the second subfunction using  Eq.36
    length_omega=length(omega);
    Lw=0;

    for i=1:length_omega
        Lw=Lw+nu*omega(i)^2-lambda(i)*omega(i);
    end
    Lw_total(iter_LR)=0;
    y=(L_low_triangular'*x_let_dynamic)'*(L_low_triangular'*x_let_dynamic); 
    
    Ly=0;%caculate Ly using Eq.34 with Min{0,zeta*sqrt(y')-nu*y'}
    Ly_total(iter_LR)=Ly;
    
    %% updata LB,UB,and epsilon using Eqs.39-41
    lowbound_optimal=Lx;
    upperbound_optimal=mu'*x_let_dynamic+zeta*sqrt(x_let_dynamic'*covSigma*x_let_dynamic);

    lowerbound=lowbound_optimal;


    upperbound=upperbound_optimal;
     


    gap_LR=(upperbound-lowerbound)/upperbound;
    %%  step 3 update the Lagrangian multipliers (lambda,nu)
    delta=0.9;%delta is a scalar chosen between 0 and 2
 
    L_low_triangular_x_let_dynamic=L_low_triangular'*x_let_dynamic;
   
  
    for i=1:num_edges
        theta_numerator(i)=delta*(upperbound-lowerbound);
        theta_denominator(i)=0;
        for z=1:num_edges
            theta_denominator(i)=theta_denominator(i)+(L_low_triangular_x_let_dynamic(z)-omega(i))^2;
        end
        theta(i)=theta_numerator(i)/theta_denominator(i);
        % using heuristic algorithm calculate the step size of each iteration k
    end
    value_omega=0;
    for i=1:length_omega
        value_omega=value_omega+omega(i)^2;
    end
    theta_nu=delta*(upperbound-lowerbound)/((value_omega-y)^2);
    nu=nu+theta_nu*(omega'*omega);% updata nu
    value_L_x_let=L_low_triangular'*x_let_dynamic;
    for i=1:num_edges
        lambda(i)=lambda(i)+theta(i)*(value_L_x_let(i)-omega(i));%updata lambda
    end
    
    iter_LR=iter_LR+1;

    
end


end










