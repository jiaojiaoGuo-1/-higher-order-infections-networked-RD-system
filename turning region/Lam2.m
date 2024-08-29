function Lambda2 = Lam2(averagek2,dE,dI,A,beta1,beta2,d,mu,gamma,averagek1)

R0 = A*beta1*averagek1/((gamma + mu)*(A + d));
R1 = A*beta2*averagek2/((A + mu)*beta1*averagek1);
delta1 = ((A + mu)*beta1*averagek1*(R1 - 1))^2 - 4*A*(A + mu)*beta1*averagek1*beta2*averagek2*(1 - R0)/R0;
I1 = ((A + mu)*beta1*averagek1*(R1 - 1) + sqrt(delta1))/(2*(A + mu)*beta2*averagek2);
E1 = ((mu - d)*I1 + d)/(A + d);
a = dE*dI;
b= -2*E1*I1*beta2*dE*averagek2 - 3*I1^2*beta2*dE*averagek2 - E1*beta1*dE*averagek1 - 2*I1*beta1*dE*averagek1 + 2*I1*beta2*dE*averagek2 + beta1*dE*averagek1 - mu*dE - dE*gamma - dI*A - dI*d;
c= 2*E1*I1*beta2*averagek2*A + 2*E1*I1*beta2*averagek2*d + I1^2*mu*beta2*averagek2 + 3*I1^2*beta2*averagek2*A + 2*I1^2*beta2*averagek2*d + E1*beta1*averagek1*A + E1*beta1*averagek1*d + I1*mu*beta1*averagek1 + 2*I1*beta1*averagek1*A + I1*beta1*averagek1*d - 2*I1*beta2*averagek2*A - 2*I1*beta2*averagek2*d - beta1*averagek1*A - beta1*averagek1*d + mu*A + mu*d + A*gamma + gamma*d;
Lambda2 = (-b+sqrt(b*b-4*a*c))/(2*a);

end