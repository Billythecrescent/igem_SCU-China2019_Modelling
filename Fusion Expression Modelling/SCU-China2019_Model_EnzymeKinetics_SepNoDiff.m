function [f] = dXdT(t, x)

k1F = 1
k1R = 0.1
k2F = 2
k2R = 0.1
k3F = 1
k3R = 0.1
k4F = 1
k4R = 0.1
E01 = 1
E02 = 1

s1=x(1);
e1=x(2);
p1=x(3);
e2=x(4);
p2=x(5);

ds1dt = -k1F*e1*s1 + k1R*(E01-e1);
de1dt = -k1F*e1*s1 + (k1R+k2F)*(E01-e1) - k2R*p1*e1
dp1dt = k2F*(E01-e1) - k2R*p1*e1 - k3F*e2*p1 + k3R*(E02-e2)
de2dt = -k3F*e2*p1 + (k3R+k4F)*(E02-e2) - k4R*p2*e2
dp2dt = k4F*(E02-e2) - k4R*p2*e2
f = [ds1dt; de1dt; dp1dt; de2dt; dp2dt];

end

tspan = 0:0.001:10;
initial = [1,1,0,1,0];

[t,x] = ode45( @dXdT, tspan, initial);
[t,x]
%plot(t,x(:,5));
%xlabel('t'); ylabel('x')