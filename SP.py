import numpy as np
import func
from func import Map
import os

def zhang2017(mymap, zeta, MaxIter, gap):
    #SP
    ## extract information a map, and perform cholesky decomposition over Sigma

    P_low_triangular = np.linalg.cholesky(mymap.sigma) #cholesky decomposition; get the lower triangular matrix

    ## initialization ; f 是path时间的均值和方差的线性组合值 f_D是(mu+P*lambda)'*X的最优值
    lmd = np.ones(mymap.n_link).reshape(-1,1)

    upper_f = 1000000 # upper_f is  upper bound for f*
    lower_fD = -1000000 # low bound for fD
    iter_SP = 1

    gap_RDG = 1 # relative inner gap
    path_rsp_SP = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s)[2]

    ## constraint generation
    paths = []
    while iter_SP < MaxIter and (upper_f-lower_fD)/upper_f >= gap:
        dynamic_mu = mymap.mu + np.dot(P_low_triangular, lmd)
        dynamic_mu = np.abs(dynamic_mu)
        paths.append(func.dijkstra(mymap.G, mymap.r_0, mymap.r_s, dynamic_mu)[2])

        # compute y and update lower bound of fD
        optimal_fD = np.dot(dynamic_mu.T, paths[-1])
    
        y = lmd.T.dot(P_low_triangular.T).dot(paths[-1])/(zeta**2) * lmd
        lower_fD = max(lower_fD, optimal_fD)
        
        # update upper bound of f
        optimal_f = np.dot(mymap.mu.T, paths[-1]) + zeta * np.sqrt(paths[-1].T.dot(mymap.sigma).dot(paths[-1]))
        if upper_f > optimal_f:
            upper_f = optimal_f
            path_rsp_SP = paths[-1]

        ## update lambda
        dual_value = np.zeros(iter_SP)
        for j in range(iter_SP):
            dual_value[j] = np.dot(mymap.mu.T, paths[j]) + zeta * np.sqrt(paths[j].T.dot(mymap.sigma).dot(paths[j]))

        approx_optimal_dual_value = np.min(dual_value)
        delta = 0.9
        eta_numerator = delta * (approx_optimal_dual_value-optimal_fD)
        eta_denominator = np.dot((y - np.dot(P_low_triangular.T, paths[-1])).T, (y - np.dot(P_low_triangular.T, paths[-1])))
        eta = eta_numerator / eta_denominator
        sub_lambda = lmd + eta*(np.dot(P_low_triangular.T, paths[-1])-y)
        value_sub_lambda = np.sqrt(np.dot(sub_lambda.T, sub_lambda))
        lmd = zeta * sub_lambda / value_sub_lambda
        # update gap_RDG

        gap_RDG = abs(upper_f-lower_fD)/upper_f
    
        iter_SP += 1

    return path_rsp_SP, iter_SP, gap_RDG

mymap = Map()
map_dir = os.getcwd() + '/Networks/'
# mymap.generate_real_map(1, map_dir)
# mymap.update_OD([1,15])
mymap.generate_simple_map('G')
print(zhang2017(mymap, 2, 200, 0.05))