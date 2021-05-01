import numpy as np
import func
import os
import time
from func import Map
from qpsolvers import solve_qp

def func_dijkstraPP(mymap, mu, edge_index):
    edge_start = np.where(mymap.M[:,edge_index]==1)[0].item()
    edge_end = np.where(mymap.M[:,edge_index]==-1)[0].item()

    path_1 = func.dijkstra(mymap.G, mymap.r_0, edge_start, mu)[2]
    path_2 = func.dijkstra(mymap.G, edge_end, mymap.r_s, mu)[2]
    
    path_dijkstraPP = (path_1 + path_2).reshape(-1)
    path_dijkstraPP[edge_index] = 1

    return path_dijkstraPP

def func_max_flow(mymap, x_t, thres=0.01):
    path_element_set = []
    mu_prime = mymap.mu.copy()

    while np.max(x_t) > thres: # if there is remaining x_rsp_quadprog left, continue;
        print("rsp", np.sum(x_t))
        # dynamic_x_rsp_quadprog = x_rsp_quadprog.copy()
        mask = x_t <= thres
        x_t[mask] = 0
        mu_prime[mask] = 100000

        _, path, path_onehot = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s, mu_prime)

        if np.min(x_t[path]) == 0:
            break

        x_t[path] -= np.min(x_t[path])
        path_element_set.append(path_onehot.reshape(-1))

        # dynamic_x_rsp_quadprog[mask] = 100

        # edge_value = np.min(dynamic_x_rsp_quadprog)
        # edge_index = np.argmin(dynamic_x_rsp_quadprog) # find the edge list with the smallest value of x_rsp_quadprog  

        # path_dijkstraPP = func_dijkstraPP(mymap, dynamic_mu, edge_index) # 找到最短路径
        # path_element_set.append(path_dijkstraPP)
        # x_rsp_quadprog -= edge_value * path_dijkstraPP # update x_rsp_quadprog

        # mask = x_rsp_quadprog < 0.01
        # x_rsp_quadprog[mask] = 0


    return path_element_set

def func_rsp_path(mymap, Q, l, t, path_rsp_GP3):
    MAX_objective_value = 10000000
    small_value = 0

    ## prepare the parameters, and solve the continous relaxation problem
    lb = np.zeros((mymap.n_link, 1)) + small_value * np.ones((mymap.n_link, 1))
    ub = np.ones((mymap.n_link, 1)) - small_value * np.ones((mymap.n_link, 1))
    Q = (Q + Q.T) / 2
    x_rsp_quadprog = solve_qp(2*Q, l.reshape(-1), G=mymap.mu.T, h=np.array(t).reshape(-1), A=mymap.M, b=mymap.b.reshape(-1), lb=lb.reshape(-1), ub=ub.reshape(-1), solver='osqp')
    # x_rsp_quadprog, exitflag = func.cvxopt_qp(2*Q, l, mymap.M, mymap.b, lb, ub) # x_rsp_quadprog is the optimal solution of the (continously) relaxed problem
    
    if x_rsp_quadprog is not None:
        # value = np.dot(mymap.mu.T, x_rsp_quadprog)
    # if  exitflag == 1 or value < t:
        dynamic_mu = mymap.mu.copy()
        ## find the set of elementary paths: path_
        path_element_set = []

        while np.sum(x_rsp_quadprog) > 0.01: # if there is remaining x_rsp_quadprog left, continue;
            print("rsp", np.sum(x_rsp_quadprog))
            dynamic_x_rsp_quadprog = x_rsp_quadprog.copy()
            mask = x_rsp_quadprog < 0.01
            dynamic_x_rsp_quadprog[mask] = 100
            dynamic_mu[mask] = 100000
            x_rsp_quadprog[mask] = 0

            edge_value = np.min(dynamic_x_rsp_quadprog)
            edge_index = np.argmin(dynamic_x_rsp_quadprog) # find the edge list with the smallest value of x_rsp_quadprog  

            path_dijkstraPP = func_dijkstraPP(mymap, dynamic_mu, edge_index) # 找到最短路径
            path_element_set.append(path_dijkstraPP)
            x_rsp_quadprog -= edge_value * path_dijkstraPP # update x_rsp_quadprog

            mask = x_rsp_quadprog < 0.01
            x_rsp_quadprog[mask] = 0
        # path_element_set = func_max_flow(mymap, x_rsp_quadprog)

        ## select the best elementary path as the solution
        path_element_array = np.array(path_element_set)
        path_mean_value = np.dot(path_element_array, mymap.mu)
        temp = np.diag(path_element_array.dot(Q.dot(path_element_array.T))).reshape(-1,1) + np.dot(path_element_array, l)
        objective_values = np.where(path_mean_value > t, MAX_objective_value, temp)

        obj_min = np.min(objective_values)
        index = np.argmin(objective_values)
        x_rsp = path_element_set[index] # 找到方差最小的路径

    else:
        x_rsp = path_rsp_GP3
        obj_min = None

    return obj_min, x_rsp.reshape(-1,1)

def func_simplepathleastvar(mymap, x_let):
    ## solve the continous relaxation problem
    lb = np.zeros(mymap.n_link)
    ub = np.ones(mymap.n_link)
    q = np.zeros(mymap.n_link)
    x_rsp_quadprog = solve_qp(2*mymap.cov, q, A=mymap.M, b=mymap.b.reshape(-1), lb=lb, ub=ub, solver='osqp')
    # x_rsp_quadprog, exitflag = func.cvxopt_qp(2*mymap.cov, p, mymap.M, mymap.b, lb, ub) # x_rsp_quadprog is the optimal solution of the (continously) relaxed problem

    if x_rsp_quadprog is not None:
        ## find the set of elementary paths: path_
        mu_dynamic = mymap.mu.copy()
        path_element_set = []

        while np.sum(x_rsp_quadprog) > 0.01: # if there is remaining x_rsp_quadprog left, continue;
            print("svar", np.sum(x_rsp_quadprog))
            x_rsp_quadprog_value = x_rsp_quadprog.copy()
            mask = x_rsp_quadprog <= 0.001
            x_rsp_quadprog_value[mask] = 100
            mu_dynamic[mask] = 100000

            edge_value = np.min(x_rsp_quadprog_value)
            edge_index = np.argmin(x_rsp_quadprog_value) # find the edge list with the smallest value of x_rsp_quadprog  

            path_dijkstraPP = func_dijkstraPP(mymap, mu_dynamic, edge_index) # 找到最短路径
            path_element_set.append(path_dijkstraPP) #save the elementary paths
            x_rsp_quadprog -= edge_value * path_dijkstraPP # update x_rsp_quadprog; decrease the original x_rsp_quadprog value
        # path_element_set = func_max_flow(mymap, x_rsp_quadprog)

        path_element_array = np.array(path_element_set)
        objective_values = np.diag(path_element_array.dot(mymap.cov.dot(path_element_array.T)))

        obj_min = np.min(objective_values)
        index = np.argmin(objective_values)
        x_rsp = path_element_set[index] # 找到方差最小的路径

    else:
        x_rsp = x_let
        obj_min = None

    return obj_min, x_rsp.reshape(-1,1)


def GP3(mymap, zeta, gap):
    ## calculate t_min and t_max (lower bound and upper bound)
    # step 1: calculate the LET path: x_let;
    x_let = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s)[2]

    # step 2: calculate the smallest variance path: x_var, and the variance
    # value (var_min)
    x_var = func_simplepathleastvar(mymap, x_let)[1]

    # step 3: calculate t_min and t_max
    t_min = np.dot(mymap.mu.T, x_let) + zeta * np.sqrt(x_var.T.dot(mymap.cov).dot(x_var)) # t_min is the sum of the two min values
    t_max = np.dot(mymap.mu.T, x_let) + zeta * np.sqrt(x_let.T.dot(mymap.cov).dot(x_let)) # t_max is one of the feasible path's value which the LET path

    gap_GP3 = abs(t_max-t_min)/t_max

    ## initialize t=(t_min+t_max)/2, and calculate matrix H, lambda_min and convert to Q and l;
    ## initialize x_rsp as x_let

    H = zeta**2 * mymap.cov - np.dot(mymap.mu, mymap.mu.T)

    lmd = np.linalg.eig(H)[0]
    lmd_min = np.min(lmd)
    Q = H - (lmd_min-1) * np.eye(mymap.n_link)
    #l = 2 * t * mymap.mu + (lmd_min-1) * np.ones((mymap.n_link, 1)) %no need to initialize l
    #here, because we will re-initialize l in the while loop

    path_rsp_GP3 = x_let
    iter_GP3 = 1

    ## loop and update t until convergence. if no feasible path, set t_min = t; if feasible path, set t_max = t;
    while gap_GP3 > gap: 
        print("gp3", gap_GP3)
        t = (t_min + t_max) / 2

        # update l, 
        l = 2 * t * mymap.mu + (lmd_min-1) * np.ones((mymap.n_link, 1))
        # find the optimal path for the constrained quadratic program
        x_rsp_propose = func_rsp_path(mymap, Q, l, t, path_rsp_GP3)[1]
        rsp_min = np.dot(mymap.mu.T, x_rsp_propose) + zeta * np.sqrt(x_rsp_propose.T.dot(mymap.cov).dot(x_rsp_propose))

        if rsp_min < t_max:
            t_max = rsp_min
            path_rsp_GP3 = x_rsp_propose

        if rsp_min < t_min:
            t_min = rsp_min
        elif rsp_min < t:
            t_max = rsp_min
            path_rsp_GP3 = x_rsp_propose
        else :
            t_min = t
        
        gap_GP3 = abs(t_max-t_min)/t_max

        iter_GP3 += 1
    
    return path_rsp_GP3, iter_GP3, gap_GP3.item()

mymap = Map()
map_dir = os.getcwd() + '/Networks/'
mymap.generate_real_map(4, map_dir)
mymap.update_OD([15, 750])
ret = GP3(mymap, 10, 0.1)
print(ret)