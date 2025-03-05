%% 1 uzd
syms x  
y(x)= (x^.2).*sin(pi - 3*x)^3 + tan(2*x+sqrt(1+x.^2)^3);

%% 2 uzd

clc, clearvars, format compact

y(x) = @(x)log((x.^5) - sin(1./x));
x_vert = 1:0.5:7; y_vert = y(x_vert);

%% 3 uzd

clc, clearvars, format compact, close all

y = @(x)5*(exp(x)).*asin(x - sqrt(x.^3).^4);
fplot(y,[0 1]);

%% 4 uzd
clc, clearvars, format compact, close all

y1 = @(x)cos(4-3*x).^3;
y2 = @(x)4 - cos(3*x).^3;
x = -pi:pi;
plot(x,y1(x),':b',x,y2(x),'-k','LineWidth',2)
legend('y1','y2'), grid on


