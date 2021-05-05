# GP3
# source code of GP3
This file shows the Python and matlab code of SP, SG, LR and GP3 algorithm proposed in the paper "Adaptive Reliable Path Planning in Gaussian Process Regulated Environment".
# File of maltlab
The maltlab folder contains the maltlab code for the four algorithms,and because the MATLAB code is too large to upload completely, only the sioux Falls road network is displayed in this folder.
# Required version of Matlab
Matlab 2018a
# Description
 - covarianceMatrix.m：method of generating covariance matrix
 - func_GP3.m and func_GP3_accelerate.m:two implementation methods of GP3
 - func_dijkstra.m and func_dijkstraPP:two implementation methods of Dijkstra
 - func_optimalpath.m:optimal path generation method
 - func_rsp_SP.m,func_rsp_zwl and func_rsp_zyl:method of SP, LR and SG
 - func_sioux_Amap.m:the network of Sioux Falls
 - main_sioux_network.m:Sample code for testing GP3 and benchmarks on Sioux Falls network

Please contact the author (guohongliang_uestc@163.com) to fetch the data for the map information, travel time statistics, as those files are too big to be uploaded.
