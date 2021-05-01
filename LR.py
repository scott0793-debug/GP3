import numpy as np
import func
from func import Map

def zeng2015(mymap, zeta, MaxIter, gap):
    ## extract information a map, and perform cholesky decomposition over Sigma

    L_low_triangular = np.linalg.cholesky(mymap.sigma) # cholesky decomposition; get the lower triangular matrix

    ## initialization, calculate lowerbound and upperbound of the problem
    lmd = np.random.rand(mymap.n_link, 1) #randomly generate the initial values for lamda,v (0,1)
    nu = np.random.rand()
    iter_LR = 1 # initialize iteration number
    x_let = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s)[2]

    y_pri = np.dot(np.dot(L_low_triangular.T, x_let).T, np.dot(L_low_triangular.T, x_let)) #the travel time variance of the shortest path
    # y' is the travel time variance of the shortest distance path

    upperbound = np.dot(mymap.mu.T, x_let) + zeta * np.sqrt(x_let.T.dot(mymap.sigma).dot(x_let)) #calculate UB using Eq.40
    lowerbound = -100
    gap_LR = (upperbound-lowerbound)/upperbound
    L_high_triangular = L_low_triangular.T

    ## step2 slove the decomposed dual problems
    while iter_LR < MaxIter and abs(gap_LR) > gap:
        ## solve a transformed shortest path problem

        lmd_L = np.dot(lmd.T, L_high_triangular).reshape(-1)
        nan_position = np.isnan(lmd_L)
        for i in range(mymap.n_link):
            if nan_position[i] == 1:
                lmd_L[i] = 0
                lmd[i] = 0

        for i in range(mymap.n_link):
            if lmd_L[i] < 0:
                lmd_L[i] = -lmd_L[i]

        dynamic_mu = mymap.mu + lmd_L.reshape(-1,1)

        path_rsp_LR = func.dijkstra(mymap.G, mymap.r_0, mymap.r_s, dynamic_mu)[2]
        
        ## calculate Lx,Ly,Lw
        
        Lx = np.dot(dynamic_mu.T, path_rsp_LR)
    
        omega = lmd/(2*nu) # calculate omega ,the second subfunction using  Eq.36
        Lw = 0

        for i in range(mymap.n_link):
            Lw = Lw + nu * omega[i]**2 - lmd[i]*omega[i]

        y = np.dot(np.dot(L_low_triangular.T, path_rsp_LR).T, np.dot(L_low_triangular.T, path_rsp_LR))
        
        Ly = 0 # caculate Ly using Eq.34 with Min{0,zeta*sqrt(y')-nu*y'}
        
        ## updata LB,UB,and epsilon using Eqs.39-41
        lowbound_optimal = Lx
        upperbound_optimal = np.dot(mymap.mu.T, path_rsp_LR) + zeta * np.sqrt(path_rsp_LR.T.dot(mymap.sigma).dot(path_rsp_LR))

        lowerbound = lowbound_optimal
        upperbound = upperbound_optimal

        gap_LR = (upperbound-lowerbound)/upperbound

        ##  step 3 update the Lagrangian multipliers (lambda,nu)
        delta = 0.9 #delta is a scalar chosen between 0 and 2
    
        L_x_let_dynamic = np.dot(L_low_triangular.T, path_rsp_LR)
    
        theta_numerator = np.zeros(mymap.n_link)
        theta_denominator = np.zeros(mymap.n_link)
        theta = np.zeros(mymap.n_link)
        for i in range(mymap.n_link):
            theta_numerator[i] = delta * (upperbound-lowerbound)
            theta_denominator[i] = 0

            for z in range(mymap.n_link):
                theta_denominator[i] = theta_denominator[i] + (L_x_let_dynamic[z]-omega[i])**2
            
            theta[i] = theta_numerator[i] / theta_denominator[i]
            # using heuristic algorithm calculate the step size of each iteration k

        value_omega = np.dot(omega.T, omega)

        theta_nu = delta*(upperbound-lowerbound)/((value_omega-y)**2)
        nu = nu + theta_nu * value_omega # updata nu
        value_L_x_let = np.dot(L_low_triangular.T, path_rsp_LR)

        for i in range(mymap.n_link):
            lmd[i] = lmd[i] + theta[i] * (value_L_x_let[i]-omega[i]) #updata lambda
        
        iter_LR += 1
        
    return path_rsp_LR,iter_LR,gap_LR

mymap = Map()
mymap.generate_simple_map('G')
print(zeng2015(mymap, 2, 200, 0.05))

