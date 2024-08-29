%参数
load('result_m.mat','result_m');
dE= 10;
dI=0.1;
load('L_ER.mat','L');
% 计算特征值和特征向量
[V, D1] = eig(L);
disp(D1);%D1中存放1000个拉普拉斯特征值

phi = 0.16;
beta1 = 0.3; 
beta2 = 0.7;
rho   = 0.25;
alpha = 0.7;
gamma = 0.1;
k1 = 7;
k2 = 15;
% 画图
figure(1)
for i=1:1000
Lambda=D1(i,i);
R1 = phi*beta2*k2/((phi + alpha)*beta1*k1);
R0 = phi*beta1*k1/((gamma + alpha)*(phi + rho));
delta1 = ((phi + alpha)*beta1*k1*(R1 - 1))^2 - 4*phi*(phi + alpha)*beta1*k1*beta2*k2*(1 - R0)/R0;
I1 = ((phi + alpha)*beta1*k1*(R1 - 1) + sqrt(delta1))/(2*(phi + alpha)*beta2*k2);
E1 = ((alpha - rho)*I1 + rho)/(phi + rho);
detk = (Lambda*dE - phi - rho)*(beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(-beta1*k1 + beta2*k2)*I1 - 3*beta2*k2*I1*I1 + dI*Lambda) - (alpha - rho)*(-I1*I1*beta2*k2 - I1*beta1*k1);
trk = Lambda*dE + Lambda*dI -(phi + rho) + (beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(beta2*k2-beta1*k1)*I1 - 3*beta2*k2*I1*I1);
x= (trk + sqrt(trk*trk - 4*detk))/2;
plot(Lambda,x,'g.');
hold on
end

k2 = 17;
% 
figure(1)
for i=1:1000
Lambda=D1(i,i);
R1 = phi*beta2*k2/((phi + alpha)*beta1*k1);
R0 = phi*beta1*k1/((gamma + alpha)*(phi + rho));
delta1 = ((phi + alpha)*beta1*k1*(R1 - 1))^2 - 4*phi*(phi + alpha)*beta1*k1*beta2*k2*(1 - R0)/R0;
I1 = ((phi + alpha)*beta1*k1*(R1 - 1) + sqrt(delta1))/(2*(phi + alpha)*beta2*k2);
E1 = ((alpha - rho)*I1 + rho)/(phi + rho);
detk = (Lambda*dE - phi - rho)*(beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(-beta1*k1 + beta2*k2)*I1 - 3*beta2*k2*I1*I1 + dI*Lambda) - (alpha - rho)*(-I1*I1*beta2*k2 - I1*beta1*k1);
trk = Lambda*dE + Lambda*dI -(phi + rho) + (beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(beta2*k2-beta1*k1)*I1 - 3*beta2*k2*I1*I1);
x= (trk + sqrt(trk*trk - 4*detk))/2;
plot(Lambda,x,'b.');
hold on
end
k2 = 19;
% 
figure(1)
for i=1:1000
Lambda=D1(i,i);
R1 = phi*beta2*k2/((phi + alpha)*beta1*k1);
R0 = phi*beta1*k1/((gamma + alpha)*(phi + rho));
delta1 = ((phi + alpha)*beta1*k1*(R1 - 1))^2 - 4*phi*(phi + alpha)*beta1*k1*beta2*k2*(1 - R0)/R0;
I1 = ((phi + alpha)*beta1*k1*(R1 - 1) + sqrt(delta1))/(2*(phi + alpha)*beta2*k2);
E1 = ((alpha - rho)*I1 + rho)/(phi + rho);
detk = (Lambda*dE - phi - rho)*(beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(-beta1*k1 + beta2*k2)*I1 - 3*beta2*k2*I1*I1 + dI*Lambda) - (alpha - rho)*(-I1*I1*beta2*k2 - I1*beta1*k1);
trk = Lambda*dE + Lambda*dI -(phi + rho) + (beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta2*k2*E1*I1 + 2*(beta2*k2-beta1*k1)*I1 - 3*beta2*k2*I1*I1);
x= (trk + sqrt(trk*trk - 4*detk))/2;
plot(Lambda,x,'k.');
hold on
end


for i=1:1000
Lambda=D1(i,i);
I1 = (beta1*k1*phi-(phi + rho)*(gamma+alpha))/(beta1*k1*(alpha+phi));
E1 = ((alpha - rho)*I1 + rho)/(phi + rho);
detk = (Lambda*dE - phi - rho)*(beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta1*k1*I1+Lambda*dI) - (alpha - rho)*(-beta1*k1*I1);
trk = (Lambda*dE - phi - rho)+(beta1*k1 - gamma - alpha - beta1*k1*E1 - 2*beta1*k1*I1+Lambda*dI);
x= (trk + sqrt(trk*trk - 4*detk))/2;
plot(Lambda,x,'m.');
hold on
end

line([min(-9) max(0)], [0 0], 'Color', 'r', 'LineStyle', '--');
dataFolder='/fig';
filename = [dataFolder,'sesan'];
saveas(gcf,filename,'fig');
print(gcf,filename,'-djpeg','-r600');

saveas(figure(1), 'sesan.eps', 'eps')