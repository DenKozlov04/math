%2 uzdevums
%atan(2*3.2-1).^3 + asin (sqrt(0.5)).^3
% definesana skaitlisko funkciju y=@(x)x^(1/3)
% y_vert_precizi=y(9)
% y_vert_precizi=y(sym(9))
%y=@(x)x.^2
%A=[1,2;3,4];
%y(A)

%3 uzdevums
%syms x;
%y(x)=(atan((sqrt(x))^3)/x^2)+((log(cot(5*x)))/(3*x*sqrt(x)));
%x_0=y(0.5)
%x_0 = double(x_0)

% 4 uzdevums
%f(x)=sqrt((cos(x)^2)+(sqrt(x^3))^4)
%x0 = 1:2:21
%vert(double(
%y0=y(x0);
%f_han = function_handle with value:@(x)cot(log(x.*2.0+1.0).^2)+log(exp(x+2.0)+asin(x.*2.0+1.0)).^3./5.0
%f_ver = -2.09727802


% 1 uzdevums , Grafiki
%syms x;
%izt=2*(atan(2*x+3)^2);
%fplot(izt,[-2*pi pi]);

%y(x)=izt;

% 2 uzdevums , Grafiki
%clc, clearvars, format compact, close all
%xt = @(t)2*cos(3*t); yt = @(t)3*sin(2*t);
%fplot(xt,yt,[-pi,pi])
%xt = matlabFunction(2*cos(3*t));
%yt = matlabFunction(3*sin(2*t));

% 3 uzdevums , Grafiki
%clc, clearvars, format compact, close all
%syms x;
%yx = @(x)3*cos(2*x.^2); 
%zx = @(x)3*(cos(2*x))^2;

%fplot(yx,[-4,4]),hold on 
%fplot(zx,[-4,4]),hold off

% 4 uzdevums , Grafiki//%% izpilda visu sekciju
clc, clearvars, format compact, close all
syms x;
yx = @(x)cos(2^x);
fplot(yx,[-2*pi,2*pi],'r-x','LineWidth', 2 ) 
title('Funkcijas 𝑦 = cos(2𝑥) grafiks')

% 5 uzdevums , Grafiki// izpilda visu sekciju

clc, clearvars, format compact, close all
syms x;
yx = @(x)(2*acot((x.^3)+3));
zx = @(x)(2*(acot(x))^3+3);


fplot(yx,[-10,10],':b','LineWidth', 2.5),hold on %2 dažadas funkcijas, tapēc 'hold on' un 'hold off'
fplot(zx,[-10,10],':r','LineWidth', 2),hold off
title('Divu funkciju grafiki')
legend('y = 2arcctg(x^3+3)','y = 2arcctg(x)^3+3')
grid on 

%% 5 uzdevums , Grafiki// izpilda visu sekciju
clc, clearvars, format compact, close all
syms t1;
syms t2;

xt1 = @(t1)(2*cos(t1)^2+cos(t));
yt1 = @(t1)(sin(2*t1)+sin(t));

xt2 = @(t2)((-cos(t)-(2*cos(0,5*t)));
yt2 = @(t2)();

fplot(yx,[-10,10],':b','LineWidth', 2.5),hold on %2 dažadas funkcijas, tapēc 'hold on' un 'hold off'
fplot(zx,[-10,10],':r','LineWidth', 2),hold off
title('Divu funkciju grafiki')
legend('y = 2arcctg(x^3+3)','y = 2arcctg(x)^3+3')
grid on 



