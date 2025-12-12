%6.1--------------------------------------------------------------------------------------------------------
% 1) Izmantojot hordu metodi, atrast vienādojuma
% 0.1x^5 + x^4 - 3x = arctg3x - 2cosX
% tuvinato sakni péc 4 iterācijām, ja sakuma tuvinājums [1 2].
% Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums1
clc; clear; format compact

f = @(x) 0.1*x.^5+x.^4-3*x-atan(3*x)+2*cos(x);
a = 1; 
b = 2;

for k = 1:4
    x = b-f(b)*(b-a)/(f(b)-f(a));
    a = b;
    b = x;
end

fprintf("%.4f\n",x);

% 2) Atrisināt vienādojumu 0.2x^3 – 5x + 5 = e^(0,3x) ar Ņūtona metodi.
% Atrast vislielāko vienādojuma sakni (precizitāte 10^(-4)).
% Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums2
clc; clear; format compact

f  = @(x) 0.2*x.^3-5*x+5-exp(0.3*x);
df = @(x) 0.6*x.^2-5+0.3*exp(0.3*x);

x = 5;
eps = 1e-4;

while true
    x_new = x-f(x)/df(x);
    if abs(x_new-x)<eps,break;
    end
    x = x_new;
end

fprintf("%.4f\n",x_new);

%3) Atrisināt sistēmu
% 	общая скобка: sin(0.5*x2 - 3) + 2.9*x1 = -2
% 	общая скобка: 0.75*x2 + cos(2.2*x1 - 2.9) + 2 = 0
% ar Ņūtona metodi apgabalā −3 ≤ x1 ≤ 2,−5 ≤ x2 ≤ −1 (precizitāte 10^-7 ).
% Aprēķināt x2 .Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums3
clc; clear; format compact

f1 = @(x1,x2) sin(0.5*x2-3)+2.9*x1+2;
f2 = @(x1,x2) 0.75*x2+cos(2.2*x1-2.9)+2;

J = @(x1,x2) [ 2.9,0.5*cos(0.5*x2 - 3);-2.2*sin(2.2*x1 - 2.9),0.75 ];
x1 = 0;  
x2 = -3;  
eps = 1e-7;

while true
    F = [-f1(x1,x2);-f2(x1,x2)];
    d = J(x1,x2)\F;
    x1_new = x1+d(1);
    x2_new = x2+d(2);
    if norm([x1_new-x1,x2_new-x2])<eps,break;
    end
    x1 = x1_new;
    x2 = x2_new;
end

fprintf("%.4f\n",x2_new);

%4) Izmantojot hordu metodi, atrast minimalo pozitivo skaitli a,
% kas apmierina vienadojumu ∫3 сверху,1 снизу e^(-0.3x)sqrt(4+ax^3)dx = 4
% (precizitate ε = 0.0001). Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums4
clc; clear; format compact

I = @(a) integral(@(x) exp(-0.3*x).*sqrt(4+a*x.^3),1,3);
F = @(a) I(a)-4;
a1 = 0.1; 
a2 = 2;
eps = 1e-4;

while abs(a2-a1)>eps
    a_new = a2-F(a2)*(a2-a1)/(F(a2)-F(a1));
    a1 = a2;
    a2 = a_new;
end

fprintf("%.4f\n",a2);


%5) Dota vienādojumu sistēma
%общая скобка:2x^5 + 3y^5 = -75
%общая скобка:(x^2)/5 + (y^2)/5 = 1
%Uzzīmēt funkciju f1(x) un f2(x) grafikus (vienā zīmēšanas logā),
%kas atbilst vienādojumu sistēmai
%общая скобка:f1(x) = 0
%общая скобка:f2(x) = 0

%%uzdevums5
clc; clear; close all; format compact

f1 = @(x,y) 2*x.^5+3*y.^5+75;
f2 = @(x,y) x.^2/5+y.^2/5-1;

figure; hold on; grid on;

fimplicit(@(x,y) f1(x,y),[-3,3,-3,3],'LineWidth',1.6);
fimplicit(@(x,y) f2(x,y),[-3,3,-3,3],'r','LineWidth',1.6);

%6.2-----------------------------------------------------------------------------------------------------------

% 1) Atrisinat vienadojumu 0.2x^3 — 5x + 5 = e^0,3x ar Ņūtona metodi.
% Atrast vismazāko vienādojuma sakni (precizitāte 10^-4).
% Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums1
clc; clear; format compact

f  = @(x) 0.2*x.^3-5*x+5-exp(0.3*x);
df = @(x) 0.6*x.^2-5+0.3*exp(0.3*x);

x = -5;
eps = 1e-4;

while true
    x_new = x-f(x)/df(x);
    if abs(x_new-x)<eps, break; end
    x = x_new;
end

fprintf("%.4f\n",x_new);

% 2) Noteikt vienādojumu sistēmas
% общая скобка: 2x^5 + 3y^5 = -75
% общая скобка: (x^2)/5 + (y^2)/5 = 1
% reālo sakņu skaitu.

%%uzdevums2
clc; clear; format compact

f1 = @(x,y) 2*x.^5+3*y.^5+75;
f2 = @(x,y) x.^2/5+y.^2/5-1;
x_max = sqrt(5);  
x_bound = x_max*1.1; 
xv = linspace(-x_bound,x_bound,400);
count = 0;

for i = 1:length(xv)
    x = xv(i);
    y1 = fzero(@(yy) f1(x,yy),0);
    if abs(f2(x,y1))<1e-3
        count = count+1;
    end
end
fprintf("%d\n",count);

% 3)Izmantojot Ņūtona metodi, atrast vienādojumu sistēmas

% общая скобка: 3y^2 - 2x^2 = 16
% общая скобка: x^3 + y^3 - 4xy = 2
% tuvinātas saknes normu ||z^(2)||_2 pēc 2 iterācijām,
% ja z^(n) = (x^(n)) un sākuma tuvinājums ir x^(0) = 2, y^(0) = -2.
%            (y^(n)) 
% Atbildi dod ar četriem cipariem aiz komata.

%%uzdevums3
clc; clear; format compact

f1 = @(x,y) 3*y.^2-2*x.^2-16;
f2 = @(x,y) x.^3+y.^3-4*x.*y-2;
J = @(x,y) [ -4*x,6*y;3*x.^2-4*y,3*y.^2-4*x ];

x = 2; 
y = -2;

for k = 1:2
    F = [-f1(x,y);-f2(x,y)];
    d = J(x,y)\F;
    x = x+d(1);
    y = y+d(2);
end

n = norm([x y],2);
fprintf("%.4f\n",n);

% 4) Noteikt vienādojuma sin 3x = (x+2)/4 - 3 reālo sakņu skaitu.

%%uzdevums4
clc; clear; format compact

f = @(x) sin(3*x)-(x+2)/4+3;
step1_1 = -1+3; 
step1_2 = step1_1*4;  
x_min = step1_2-2;  
step2_1 = 1+3;  
step2_2 = step2_1*4; 
x_max = step2_2-2;

x = linspace(x_min, x_max, 20000);
fx = f(x);
count = 0;

for k = 1:length(x)-1
    if fx(k)*fx(k+1) < 0
        count = count + 1;
    end
end

fprintf('%d\n', count);

% 5) Dots vienadojums 2ln(x + 3) + x^2 = -3x + 5.
% Uzzimet funkcijas f(x) grafiku, kas atbilst vienadojumam f(x) = 0.

%%uzdevums5
clc; clear; close all; format compact

f = @(x) 2*log(x+3)+x.^2+3*x-5;
fplot(f,[-2.9,5],'LineWidth',1.6);grid on;

7.1--------------------------------------------------------------------------------------------------------------
%+++++++++
% 1) Atrisinat Košī problēmu (1 + 3x)y' + 2xy + (1 + 3e^x)y'' = 4 + x^2，
% y(0.22) = 0.38 , y'(0.22) = 1.55, intervālā [ 0.22, 6.77].
% Uzzimét funkcijas y'(x) grafiku.

%1uzdevums
clc; clearvars; close all;
format compact;

x0 = 0.22;
xf = 6.77;
y0 = 0.38;
yp0 = 1.55;

odefun = @(x,z) [z(2);(4+x.^2-(1+3*x).*z(2)-2*x.*z(1))./(1+3*exp(x))];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(odefun,[x0,xf],[y0;yp0],opts);

xx = linspace(x0,xf,600);
zz = deval(sol,xx);
yprime = zz(2,:);

figure;
plot(xx,yprime,'LineWidth',1.2);
grid on;

%+++++++++
% 2) Atrisināt Košī problēmu y' - x^2siny = -y^2x,
% y(0) = 1 intervālā [0, 2.5].
% Aprēķināt y(1.77) ar precizitāti 0.0001.
% Atbildi dod ar četriem cipariem aiz komata.


%2uzdevums
clc; clearvars; close all;
format long;

x0 = 0;
xf = 2.5;
x_query = 1.77;
y0 = 1;

odefun = @(x,y)x.^2.*sin(y)-x.*y.^2;
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],y0,opts);

y_at = deval(sol,x_query);
fprintf('y(%.2f) = %.4f',x_query,y_at);

%+++++++++
% 3)Atrisināt Košī problēmu
% dy1/dx = 1/2 * y1^2 + x * y2
% dy2/dx = -1/2 * y2^2 + y1 * y2
% y1(0.2) = 1.25, y2(0.2) = 1.33 intervālā (0.2, 0.9).
% Uzzīmēt funkcijas y2(x) grafiku.

%3uzdevums
clc; clearvars; close all;
format compact;

x0 = 0.2;
xf = 0.9;
y1_0 = 1.25;
y2_0 = 1.33;

odefun = @(x,y)[0.5*y(1).^2+x.*y(2); -0.5*y(2).^2+y(1).*y(2)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(odefun,[x0,xf],[y1_0;y2_0],opts);

xx = linspace(x0,xf,400);
YY = deval(sol,xx);
y2 = YY(2,:);

figure;
plot(xx, y2,'LineWidth',1.2);
grid on;

%-------------------------
%4) Atrisināt Košī problēmu
%((2sinx)/x + arctg(3 + x))y'' + sqrt(5 + 3x^3)y'+ y''' + ln(cosx + x^2)y = 5x^3 + 2,
%y(2.23) = 3.55, y'(2.23) = 4.38, y''(2.23) = 7.46 intervālā [2.23, 7.42].
%Aprēķināt y''(3.41) ar precizitāti ε = 0.0001.
%Atbildi dod ar četriem cipariem aiz komata.

%4uzdevums
clc; clearvars; close all;
format long;

x0 = 2.23; 
xf = 7.42;
x_query = 3.41;

z1_0 = 3.55;   
z2_0 = 4.38;   
z3_0 = 7.46;   

odefun = @(x,z) [z(2);z(3);5*x.^3+2-((2*sin(x)./x)+atan(3+x)).*z(3)-sqrt(5+3*x.^3).*z(2)-log(cos(x)+x.^2).*z(1)];
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],[z1_0;z2_0;z3_0],opts);

z_at = deval(sol,x_query);
ypp_at = z_at(3);
fprintf('%.4f\n',x_query,ypp_at);

%-------------------------
% 5) Atrisināt Košī problēmu  

% x'(t) = 2x(t) - 3t * y(t)  
% y'(t) = 1/5 y(t) + 4x(t)  

% x(5.5) = 1, y(5.5) = 1 intervālā (5.5  8.3).  
% Atrast x(6.75). Atbildi dod ar četriem cipariem aiz komata.

%5uzdevums
clc; clearvars; close all;
format long;

t0 = 5.5;
tf = 8.3;
t_query = 6.75;
x0 = 1;
y0 = 1;

odefun = @(t,u)[2*u(1)-3*t.*u(2);4*u(1)+(1/5)*u(2)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
sol = ode45(odefun,[t0,tf],[x0;y0],opts);

u_at = deval(sol,t_query);
x_at = u_at(1);
fprintf('%.4f\n',t_query,x_at);

7.2--------------------------------------------------------------------------------------------------------------
%+++++++
% 1) Atrisināt Koši problēmu
% ((2sinx)/x + arctg(3 + x)) y'' + sqrt(5 + 3x^3)y' + y''' + ln(cos x + x^2)y = 5x^3 + 2,
% y(2.23) = 3.55, y'(2.23) = 4.38, y''(2.23) = 7.46 intervālā [2.23, 7.42].

% Uzzimēt funkcijas y'(x) grafiku.

%1uzdevums
clc; clearvars; close all;
format compact;

x0 = 2.23;
xf = 7.42;
IC = [3.55;4.38;7.46];

odefun = @(x,z) [z(2);z(3);5*x.^3 + 2-((2*sin(x)./x)+atan(3+x)).* z(3)-sqrt(5+3*x.^3).* z(2)-log(cos(x)+x.^2).* z(1)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],IC,opts);

xx = linspace(x0,xf,600);
zz = deval(sol,xx);   
yprime = zz(2,:);

figure;
plot(xx, yprime, 'LineWidth', 1.2);
grid on; 

%++++++
% 2) Atrisināt Košī problēmu (1 + sinx)y'' + (1 + x^2)y + 2y' = 8x^2,
% y(0.33) = 0.52, y'(0.33) = 2.33, intervalā [0.33, 4.55].

% Aprēķināt y'(1.99) ar precizitāti 0.0001.
% Atbildi dod ar četriem cipariem aiz komata.

%2uzdevums
clc; clearvars; close all;
format long;

x0 = 0.33;
xf = 4.55;
xq = 1.99;
IC = [0.52;2.33];

odefun = @(x,z) [ z(2);(8*x.^2-(1+x.^2).*z(1)-2.*z(2))./(1+sin(x))];
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],IC,opts);
zq = deval(sol,xq);

yprime_q = zq(2);
fprintf('%.4f\n',xq,yprime_q);


%---------
% Atrisināt Košī problēmu  

% y'''+( 5/(2x) + ln(sinx + sqrt(x^3)))y''+ sqrt(1 + x^2)y' + ln(5 + x^3)y = sqrt(x + 4),  
% y(3.55) = -1.25, y'(3.55) = 2.05, y''(3.55) = 3.42 intervalā [3.55, 7.24].  

% Aprēķināt y(5.32) ar precizitāti ε = 0.0001.  
% Atbildi dod ar četriem cipariem aiz komata.

%3uzdevums
clc; clearvars; close all;
format long;

x0 = 3.55;
xf = 7.24;
xq = 5.32;
IC = [-1.25;2.05;3.42];

odefun = @(x,z) [z(2);z(3);sqrt(x + 4)-((5./(2.*x))+log(sin(x)+x.^(3/2))).* z(3)-sqrt(1+x.^2).*z(2)-log(5+x.^3).*z(1)];

xx_check = linspace(x0,xf,1000);
arg1 = sin(xx_check)+xx_check.^(3/2);
arg2 = 5+xx_check.^3;
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

sol = ode45(odefun,[x0,xf],IC,opts);
zq = deval(sol,xq);
y_q = zq(1);
fprintf('%.4f\n', xq, y_q);

%++++++
% 4) Atrisināt Koši problēmu y' - x^2siny = -y^2x,  
% y(0) = 1 intervālā [0, 2.5].  
% Uzzīmēt funkcijas y(x) grafiku.

%4uzdevuns
clc; clearvars; close all;
format compact;

x0 = 0;
xf = 2.5;
IC = 1;

odefun = @(x,y) x.^2.*sin(y)-x.*y.^2;
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],IC,opts);

xx = linspace(x0,xf,500);
yy = deval(sol,xx);

figure;
plot(xx,yy,'LineWidth',1.2);
grid on;

%++++++++
% 5) Atrisināt Košī problēmu  

% dy1/dx = 1/3 * y1 - 0.02y1y2  
% dy2/dx = -y1 + 0.02y1y2  

% y1(1.8) = 1.25, y2(1.8) = 0.88 intervālā (1.8, 8.2).  
% Atrast y2(4.42). Atbildi dod ar četriem cipariem aiz komata.


%5uzdevums
clc; clearvars; close all;
format long;

x0 = 1.8;
xf = 8.2;
xq = 4.42;
IC = [1.25;0.88];

odefun = @(x,y) [(1/3)*y(1)-0.02*y(1)*y(2); -y(1)+0.02*y(1)*y(2)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
sol = ode45(odefun,[x0,xf],IC,opts);
Yq = deval(sol,xq);
y2_q = Yq(2);

fprintf('%.4f\n',xq,y2_q);




