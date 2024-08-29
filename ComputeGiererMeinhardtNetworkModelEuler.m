
function []= ComputeGiererMeinhardtNetworkModelEuler(N, networkType,lab,local_Cluster_Coefficient)

% sigma 扩散比例系数; N 网络节点数.
% networkType 网络类型: 'LA', 'WS', 'ER', 'BA'.
% averageDegree 网络平均度.

global d1 d2 G E

%% 模型参数
d1 = 0.01;
d2 = 1;
G = 3;

time = 500;

%% 平衡点
U_star = E/G;
V_star = E/G;

%% 图灵失稳条件的判断
a11 = 1;
a12 = - 1;
a21 = 2*E;
a22 = - E;

J = [a11,a12; a21,a22]; % 平衡点处的Jacobi矩阵
trJ = trace(J);
detJ = det(J);

disp('=============================正平衡点===============================')
fprintf('U* = %f,  V* = %f\n', U_star, V_star);
disp('==============================Turing分支条件========================')
fprintf('tr(J) = %f,  det(J) = %f\n', trJ, detJ);
fprintf('a11d2 + a22d1 = %f,  2sqrt(d1d2det(J)) = %f\n', ...
    a11*d2 + a22*d1, 2*sqrt(d1*d2*detJ));

%% 从文件读入网络的Laplacian矩阵
L = lab;
%% 调用ode求解器计算网络上的反应扩散方程

dataFolder = '.\Data\';
if ~isfolder(dataFolder) 
    mkdir(dataFolder); 
end

filename = [dataFolder,'RandnData',num2str(N),'.mat'];
if ~exist(filename,'file')
    U0 = randn(N,1);
    V0 = randn(N,1);
    save(filename,'U0','V0');
else
    load(filename,'U0','V0');
end

% 模型初值: 平衡点的随机扰动(每次计算扰动值固定)
U0 = U_star + 0.0005*U0;
V0 = V_star + 0.0005*V0;

deltaT = 0.005;
timeGrid = (0:deltaT:time)';
t = (0:1:time)';

U = zeros(length(t), N);
V = zeros(length(t), N);
U(1,:) = U0';
V(1,:) = V0';
index = 1;

for n = 2:length(timeGrid)
    U1 = U0 + f(U0,V0)*deltaT + d1*L*U0*deltaT;
    V1 = V0 + g(U0,V0)*deltaT + d2*L*V0*deltaT;
    U0 = U1;
    V0 = V1;
    
    if (mod(timeGrid(n),1) == 0)
        index = index + 1;
        U(index,:) = U0';
        V(index,:) = V0';
        fprintf('----------时间 = %f 完成计算!----------\n', timeGrid(n));
    end
end

% 按节点聚类系数降序排列U的值
U = sortrows([local_Cluster_Coefficient, U'], 1);
U = U(:,2:end)';
V = sortrows([local_Cluster_Coefficient, V'], 1);
V = V(:,2:end)';

fig = figure;
set(fig,'visible','on');

x = 1:1:N;

% figure;
plot(x, U(end,:), 'k.', 'MarkerSize', 18);
hold on;
plot(x, U_star*ones(N,1), 'r');
xlabel('i');
ylabel('u_i');
axis([0,N,-inf,inf]);
set(gca, 'FontSize',20);
set(get(gca,'Children'),'linewidth',3.0);
set(get(gca,'XLabel'),'FontSize',25);
set(get(gca,'YLabel'),'FontSize',25);
filename = [dataFolder,networkType,'_u_1D'];
saveas(gcf,filename,'fig');
print(gcf,filename,'-djpeg','-r600');

[x,t] = meshgrid(x,t);

filename = [dataFolder,networkType,'_data','.mat'];
save(filename, 'x', 't', 'U', 'V');

% figure;
% surf(x,t,U);
% colorbar;
% colorbar('FontSize',18);
% colormap jet
% xlabel('i');
% ylabel('t');
% set(gca, 'FontSize',20);
% set(get(gca,'XLabel'),'FontSize',20);
% set(get(gca,'YLabel'),'FontSize',20);
% shading interp;
% view([0 90]);
% axis([0,N,-inf,inf]);
% filename = [dataFolder,networkType,'_u_2D'];
% saveas(gcf,filename,'fig');
% print(gcf,filename,'-djpeg','-r600');

end


%% 非线性反应项GM模型
function val = f(u,v)

val = u.^2./v - u;

end

function val = g(u,v)

global G E

val = G*u.^2 - E*v;

end


