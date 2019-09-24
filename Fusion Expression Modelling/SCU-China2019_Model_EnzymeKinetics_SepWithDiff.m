function [f] = dXdT(t, x)

k1F = 1
k1R = 0.1
k2F = 2
k2R = 0.1
k3F = 1
k3R = 0.1
k4F = 1
k4R = 0.05
E1 = 1
E2 = 1

s1=x(1);
e1=x(2);
pl=x(3);
pr=x(4);
e2=x(5);
p2=x(6);

kdiff = 10**(-2);
d = 1.118/((e1+e2)**(1/3));
%d = 1.118/2**(1/3)

ds1dt = -k1F*e1*s1 + k1R*(E1-e1);
de1dt = -k1F*e1*s1 + (k1R+k2F)*(E1-e1) - k2R*pl*e1;
dpldt = k2F*(E1-e1) - k2R*pl*e1 - (pl-pr)/exp(kdiff*d);
dprdt = (pl-pr)/exp(kdiff*d) - k3F*e2*pr + k3R*(E2-e2);
de2dt = -k3F*e2*pr + (k3R+k4F)*(E2-e2) - k4R*p2*e2;
dp2dt = k4F*(E2-e2) - k4R*p2*e2;
f = [ds1dt; de1dt; dpldt; dprdt; de2dt; dp2dt];

end

tspan = 0:0.001:40;
initial = [1,1,0,0,1,0];

[t,x] = ode45( @dXdT, tspan, initial);
[t,x]
%plot(t,x(:,6));
%xlabel('t'); ylabel('x')
