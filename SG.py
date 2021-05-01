import numpy as np
import func
from func import Map

def zhang2019(mymap, zeta, MaxIter, gap):
    #SG
    #relative duality gap threshold epsilon, MaxIter is maximum number of iteration 

    ## step 1,variance-covariance maxtrix decomposition
    eig_value, eig_vector = np.linalg.eig(mymap.cov)

    ## step2 initialization
    iter_SG = 0
    lower_bound_Z = -1000
    lmd = np.ones(mymap.n_link).reshape(-1,1)
    x_let = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s)[2]

    upper_bound_Z = np.dot(mymap.mu.T, x_let) + zeta * np.sqrt(x_let.T.dot(mymap.cov).dot(x_let)) #compute Z as the upper bound
    gap_SG = (upper_bound_Z-lower_bound_Z) / upper_bound_Z

    ## step 3 solve lagtangian telaxation problems
    while iter_SG < MaxIter and gap_SG > gap:
        dynamic_mu = mymap.mu + np.dot(eig_vector, lmd)
        dynamic_mu = np.abs(dynamic_mu)
        x_let = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s, dynamic_mu)[2]
        
        optimal_upper_bound_Z = np.dot(mymap.mu.T, x_let) + zeta * np.sqrt(x_let.T.dot(mymap.cov).dot(x_let)) #update upper bound of Z
        lagrangian_optimal = np.dot((mymap.mu + np.dot(eig_vector, lmd)).T, x_let)

        ## update upper_bound_Z and save x_let
        if upper_bound_Z >= optimal_upper_bound_Z:
            upper_bound_Z = optimal_upper_bound_Z
            x_SG = x_let
            x_SG_last = x_SG
        else:
            if iter_SG == 0:
                x_SG = x_let
                x_SG_last = x_SG
            else:
                x_SG = x_SG_last

        lagrangian_value = lagrangian_optimal

        ## update low_bound_Z
        if lagrangian_value >= lower_bound_Z:
            lower_bound_Z = lagrangian_value

        # update lagrangian multiplier
        partial_coefficient = np.dot(eig_vector.T, x_let)
        coefficient = np.abs(upper_bound_Z-lower_bound_Z) / np.dot(partial_coefficient.T, partial_coefficient)
        beta = 0.1
        lmd_new = lmd + beta * coefficient * partial_coefficient
        for i in range(mymap.n_link):
            if lmd_new[i] > np.abs(zeta * np.sqrt(eig_value[i])):
                lmd_new[i] = lmd[i]
        
        lmd = lmd_new
        iter_SG += 1
        gap_SG = (upper_bound_Z-lower_bound_Z)/upper_bound_Z

    return x_SG, iter_SG, gap_SG.item()

mymap = Map()
mymap.generate_simple_map('G')
# print(zhang2019(mymap, 2, 200, 0.05))
# import networkx as nx

# graph = nx.DiGraph()
# graph.add_node(1)
# graph.add_edge(1,2)
# graph.add_node(2)
# graph.add_node(2)

# print(graph.nodes)