%% 网络上的反应扩散方程
function dy = GiererMeinhardtNetworkModel(t,y,L)

global d1 d2

N = size(L,1);
dy = zeros(2*N,1);

for i = 1:N
    dy(i) = f(y(i),y(N+i)) + d1*L(i,:)*y(1:N,1);
    dy(N+i) = g(y(i),y(N+i)) + d2*L(i,:)*y(N+1:end,1);
end


end


%% 非线性反应项
function val = f(u,v)

val = u^2/v - u;

end

function val = g(u,v)

global G E

val = G*u^2 - E*v;

end




