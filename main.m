clc; 
clear;
global A d mu beta1 beta2 gamma averagek1 averagek2 dE dI 

%% 模型参数
averagek2 =11;%

A = 0.21;
d = 0.1;
mu = 0.9;
gamma = 0.3;
beta1 = 0.4; 
averagek1 = 5;
beta2 = 0.6;


dE=10;
dI=0.1;
N = 200;

load('L_ER3.mat','L')
ComputeHONetworkModelEuler(N,L);



