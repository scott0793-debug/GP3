function [optimal_path,optimal_value]=func_optimalpath(mu,covSigma,zeta)
all_path_value=zeros(50,1);
optimal_path=zeros(50,1);
for i=1:50
    path=zeros(50,1);
    path(i)=1;
    all_path_value(i)=path'*mu+zeta*sqrt(path'*covSigma*path);
end
[optimal_value,index]=min(all_path_value);
optimal_path(index)=1;
end