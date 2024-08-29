d1= 10;
d2=0.1;

A = 0.21;
beta1 = 0.4;
gamma = 0.3;
d = 0.1;
mu = 0.9;
averagek1 = 5;
beta2 = 0.8;

p=zeros(110,1);
q=zeros(110,1);
n=0;
%% 离散
for k2=9:0.02:28
Lambda1=Lam1(averagek2,dE,dI,A,beta1,beta2,d,mu,gamma,averagek1);
Lambda2=Lam2(averagek2,dE,dI,A,beta1,beta2,d,mu,gamma,averagek1);

n=n+1;
p(n,1)=Lambda1;
q(n,1)=Lambda2;

x_points = k2;
y_points = Lambda1;%下限
plot(x_points, y_points, 'b.','LineWidth',2);
hold on
x_points = k2;
y_points = Lambda2;%上限
plot(x_points, y_points, 'b.','LineWidth',2);
hold on
end
