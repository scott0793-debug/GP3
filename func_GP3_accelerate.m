function x=func_GP3_accelerate(inv_Q,lambda_min,lambda01,mu,lambda02,lambda03,omega,A,t_min,t_max,b,num_edges,num_nodes,covSigma)
on=ones(num_edges,1);
for i=1:50
    t=(t_min+t_max)/2;
    l = 2*t*mu+(lambda_min-1)*ones(num_edges,1); 
    x=-inv_Q*(l+lambda01*mu+(omega*A)'-lambda02'+lambda03');
    lambda01=lambda01+0.1*(mu'*x-t);
    object01=A*x;
     for j=1:num_nodes
         omega(j)=omega(j)+0.1*(object01(j)-b(j));      
     end
    for j=1:num_edges
         lambda03(j)=lambda03(j)+0.1*(x(j)-on(j));
    end
    for j=1:num_edges
         lambda03(j)=lambda03(j)-0.1*(x(j));
    end
    val=mu'*x+sqrt(x'*covSigma*x);
    if val>t
        t_min=t;
    else
        t_max=val;
    end

end

end