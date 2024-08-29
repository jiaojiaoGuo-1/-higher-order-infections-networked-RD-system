function val = g(E,I)
global mu beta1 beta2 gamma averagek1 averagek2  

val = (beta1*averagek1-gamma-mu).*I-beta1*averagek1.*E.*I-beta2*averagek2.*E.*I.*I+(beta2*averagek2-beta1*averagek1).*I.*I-beta2*averagek2.*I.*I.*I;
end