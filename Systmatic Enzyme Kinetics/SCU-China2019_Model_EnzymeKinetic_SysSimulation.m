function [f] = dXdT(t, x)

k1F = 1;
k1R = 0.1;
k2F = 2;
k2R = 0.1;
k3F = 1;
k3R = 0.1;
k4F = 1;
k4R = 0.1;
k5F = 1;
k5R = 0.5;
k6F = 1;
k6R = 0.5;
k7F = 3;
k7R = 0.1;
k8F = 3;
k8R = 0.1;
k9F = 10;
k9R = 0.1;
E01 = 1;
E02 = 1;
E03 = 1;
E04 = 1;

s1=x(1);
e1=x(2);
p1=x(3);
e2=x(4);
p2=x(5);
e3=x(6);
p3=x(7);
s4=x(8);
e4=x(9);
i=x(10);
ei=x(11);

%ds1dt = -k1F*e1*s1 + k1R*(E01-e1);
dsdt = 0;
de1dt = -k1F*e1*s1 + (k1R+k2F)*(E01-e1) - k2R*p1*e1;
dp1dt = k2F*(E01-e1) - k2R*p1*e1 - k3F*e2*p1 + k3R*(E02-e2);
de2dt = -k3F*e2*p1 + (k3R+k4F)*(E02-e2) - k4R*p2*e2;
dp2dt = k4F*(E02-e2) - k4R*p2*e2 - k5F*e3*p2 + k5R*(E03-e3-ei);
de3dt = -k5F*e3*p2 + (k5R+k6F)*(E03-e3-ei) - k6R*p3*e3 - k9F*i*e3 + k9R*ei;
dp3dt = k6F*(E03-e3-ei) - k6R*p3*e3;
%ds4dt = -k7F*e4*s4 + k7R*(E04-e4);
ds4dt = 0;
de4dt = -k7F*e4*s4 + (k7R+k8F)*(E04-e4) - k8R*e4*s4;
didt  = k8F*(E04-e4) - k8R*e4*s4 - k9F*i*e3 + k9R*ei;
deidt = k9F*i*e3 - k9R*ei;
f = [ds1dt; de1dt; dp1dt; de2dt; dp2dt; de3dt; dp3dt; ds4dt; de4dt; didt; deidt];

end

tspan = 0:0.001:100;
initial = [1,1,0,1,0,1,0,2,1,0,0];

[t,x] = ode45( @dXdT, tspan, initial);
[x(:,1),x(:,5),x(:,7)]
%plot(t,x(:,5));
%xlabel('t'); ylabel('x')