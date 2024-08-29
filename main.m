clc; clear;
global E

E = 9; % 9, 10.5, 12
N = 1000;

filename = 'ClusteredNetwork9_25.mat';
load(filename,'lab','local_Cluster_Coefficient');
networktype='ER9_25';
ComputeGiererMeinhardtNetworkModelEuler(N,networktype,lab,local_Cluster_Coefficient);

