function []= ComputeHONetworkModelEuler(N,L)
global A d mu beta1 beta2 gamma averagek1 averagek2 dE dI  time

time =2000;

%% 
R1 = A*beta2*averagek2/((A + mu)*beta1*averagek1);
R0 = A*beta1*averagek1/((gamma + mu)*(A + d));
delta1 = ((A + mu)*beta1*averagek1*(R1 - 1))^2 - 4*A*(A + mu)*beta1*averagek1*beta2*averagek2*(1 - R0)/R0;
I_star = ((A + mu)*beta1*averagek1*(R1-1)+sqrt(delta1))/(2*(A + mu)*beta2*averagek2);
E_star = ((mu-d)*I_star+d)/(A + d);

%% 
a11 = -(A+d);
a12 =(mu-d);
a21 = -beta1*averagek1*I_star-beta2*averagek2*I_star*I_star;
a22 = (beta1*averagek1-gamma-mu)-beta1*averagek1*E_star-2*beta2*averagek2*I_star*E_star+2*(beta2*averagek2-beta1*averagek1)*I_star-3*beta2*averagek2*I_star*I_star;

J = [a11,a12; a21,a22]; 
trJ = trace(J);
detJ = det(J);

disp('=============================平衡点===============================')
fprintf('I* = %f,  E* = %f\n', I_star, E_star);

%% 
dataFolder = '.\Data\';
filename = [dataFolder,'RandnData',num2str(N),'.mat'];
E0 = randn(N,1);%Initial value
I0 = randn(N,1);

% 
I0 = I_star + 0.0005*I0;
E0 = E_star + 0.0005*E0;

deltaT = 0.01;
timeGrid = (0:deltaT:time)';
t = (0:1:time)';

E = zeros(length(t), N);
I = zeros(length(t), N);
E(1,:) = E0';
I(1,:) = I0';
index = 1;

for n = 2:length(timeGrid)
    E1 = E0 + f(E0,I0)*deltaT + dE*L*E0*deltaT;%????涓轰粈涔堟病鏁版嵁
    I1 = I0 + g(E0,I0)*deltaT + dI*L*I0*deltaT;
    E0 = E1;
    I0 = I1;
    
    if (mod(timeGrid(n),1) == 0)
        index = index + 1;
        E(index,:) = E0';
        I(index,:) = I0';
        fprintf('----------时间 = %f 计算完成!----------\n', timeGrid(n));
    end
end


fig = figure;
set(fig,'visible','on');

x = 1:1:N;
ylim([0 0.2]);
% figure;
plot(x, I(end,:), 'b.', 'MarkerSize', 10);
hold on;
%plot(x, I_star*ones(N,1), 'r');
xlabel('i');
ylabel('I_i');
axis([0,N,0,1]);
set(gca, 'FontSize',20);
set(get(gca,'Children'),'linewidth',3.0);
set(get(gca,'XLabel'),'FontSize',25);
set(get(gca,'YLabel'),'FontSize',25);
filename = [dataFolder,'ER3beta2_06k2','_',num2str(averagek2),'_', num2str(time)];
saveas(gcf,filename,'fig');
print(gcf,filename,'-djpeg','-r600');

[x,t] = meshgrid(x,t);

file_name= [dataFolder,'ER3beta2_06k2','_',num2str(averagek2),'_', num2str(time),'_data','.mat'];%绗嚑鍒梍鏃堕棿
save(file_name, 'x', 't', 'I', 'E');

end



