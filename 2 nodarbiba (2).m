%% 1. piemērs. fplot-funkcija uzdota atklātā veidā
clc, clearvars, format compact
% a.variants -funkcija definēta (function handle)
y = @(x)cos(2*x).^2;
fplot(y,[-pi,2*pi])

%% 2.piemērs. fplot - funkcija uzdota paramet. eidā
clc, clearvars, format compact
% b.variants - funkcija definēta (symbolic function )
syms t
xt (t) = 2*cos(t)^3;
yt (t) = 2*sin(t)^3;
fplot(xt,yt,[-pi,pi])
%% 3.piemērs. fplot. Uzzīmēt funkcijas grafiku,
% izmantojot līnijas veidu un biezumu
clc, clearvars, format compact
% a.variants - funkcija definēta ( function handle )
y = @(x)cos(x.^3+3).^2;
fplot(y,[-pi,pi],'-r', 'LineWidth',3)
%% 4.piemērs. fplot . Uzzīmēt funkciju grafikus vienā
% logā, izmantojot līnijas veidu un biezumu
% b.variants - funkcija definēta( symbolic function )
clc, clearvars, format compact, close all,syms x t
f1(x) = sin(2*x)^2; f2(x) = cos(x+2);
xt(t) = cos(3*t); yt(t) = sin(2*t);
fplot(f1,[0,2*pi],'-r','LineWidth',2),hold on
fplot(f2,[0,2*pi],'--b','LineWidth',3)
fplot(xt,yt,[0,2*pi],':g','LineWidth',4), hold off
legend('f1(x)', 'f2(x)', 'f3(x)')
title('Trigonometriskās funkcijas'), grid on
%% 5. piemērs. plot. Funkcija uzdota atklātā veidā
clc, clearvars, format compact, close all
% formula ( sākumā definē vektoru x ) 
x = -3:0.01:3; % x - vērtību vektors
y = x.^2+x.^3; % y - vērtību vektors
plot(x,y,':r','LineWidth',3), grid on
%% 6.piemērs. plot- funkcija uzdota atklātā veidā
clc, clearvars, format compact, close all
% a.variants – funkcija definēta ( function handle )
%  vektoru x nav obligāti definēt sākumā 
y = @(x)x.^2.*log(x.^3+3);
x_pr = -1:0.01:6;
plot(x_pr,y(x_pr),'-g','LineWidth',2)
grid on
%% 7.piemērs. Plot. Funkcija uzdota parametriskā veidā.
clc, clearvars, format compact , close all
% b.variants – funkcija definēta ( symbolic function)
% vektoru t var definēt pēc funkcijas uzdošanas
syms t, x(t) = sqrt(2)*cos(t)^3;
y(t) = sqrt(2)*sin(t)^3;
t = -pi:0.01:pi;  
plot(x(t),y(t), ':r','LineWidth',3)
grid on, title('Parametriskā veidā dota funkcija - Astroīda')
%% 8.piemērs. Plot. Uzzīmēt dažu fun.grafikus
clc, clearvars, format compact, clear all
% a.variants –funkcija definēta (function
% handle) vektoru x var uzdot pēc funkcijām
y1 = @(x)sin(x); y2 = @(x)sin(2*x).*cos(x);
y3 = @(x)cos(x).^2;
x = -pi:0.01:2*pi;
plot(x,y1(x),'-r',x,y2(x),'--g',... 
    x,y3(x),':b','LineWidth',3)
legend('y = sinx','y = sin2xcosx','y = cos^2x')
grid on,xlabel('x-ass'), ylabel('y-ass')
title('Dažu funkciju grafiki')
% turpinājums ( katrai līnijai ir savs biezums )
figure
zim = plot(x,y1(x),'-r',x,y2(x),'--g',x,y3(x),':b')
set(zim(1),'LineWidth',2)
set(zim(2),'LineWidth',3)
set(zim(3),'LineWidth',4)
legend('y = sinx','y = sin2xcosx','y = cos^2x')

%% Uzdevumi patstāvīgai risināšanai
%% 1.uzdevums
clc, clearvars, format compact
y = @(x)2*atan(2*x+3).^2
fplot(y,[-2*pi, pi])
%% 2.uzdevums
clc, clearvars, format compact, close all
x = @(t) 2*cos(3*t)
y = @(t) 3*sin(2*t)
t = -pi:0.01:pi;
plot(x(t),y(t))
%% 3.uzdevums
clc, clearvars, format compact, close all
y = @(x) 3*cos(2*x.^2)
z = @(x) 3*cos(2*x).^2
x = -4:0.01:4;
plot(x, y(x), x, z(x))
%% 4. uzdevums
clc, clearvars, format compact, clear all
y = @(x) cos(2.^x)
fplot(y,[-2*pi, 2*pi], '-rx', 'LineWidth', 2)
title('Funkcijas 𝑦 = cos(2^x) grafiks')
xlabel('x-ass') 
ylabel('y-ass')
%% 5.uzdevums
clc, clearvars, format compact, clear all
y = @(x) 2*acot(x.^3+3)
z = @(x) 2*(acot(x)).^3+3
x = -10:0.01:10;
zim = plot(x, y(x),'--',x,z(x),':r')
set(zim(1),'LineWidth',3)
set(zim(2),'LineWidth',2)
title('Divu funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('y = 2arccrg(x^3+3)','y = 2arccrg^3(x)+3')
%% 6.uzdevums
clc, clearvars, format compact, clear all
syms t
x1 = @(t) 2*(cos(t)).^2 + cos(t)
y1 = @(t) sin(2*t) + sin(t)
x2 = @(t) -cos(t) - 2*cos(0.5*t)
y2 = @(t) -sin(t) + 2*sin(0.5*t)
t_pr = -pi:0.01:pi;
zim = plot(x1(t_pr),y1(t_pr),'g',x2(t_pr),y2(t_pr),'b')
set(zim(1),'LineWidth',3)
set(zim(2),'LineWidth',4)
title('Divu parametriski doto funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('x1(t)y1(t)','x2(t)y2(t)')
%% 7.uzdevums
clc, clearvars, format compact, close all
syms x
f1(x) = sqrt((sin(x+2*x^2)).^3 + (3*x^5)^(1/3))
f2(x) = log(x^2 + 1) + sqrt(x)
x_ver = 2:1:12;
y_ver = double(f1(x_ver));
tabula = [x_ver' y_ver'];
disp('x_mezgli y_mezgli')
disp( tabula)

hold on
plot(x_ver, y_ver, 'Ok', 'LineWidth', 3)
fplot(f2, [2, 12], 'r-', 'LineWidth', 3)
legend('y = f_1(x)', 'y = f_2(x)')
grid on,xlabel('x-ass'), ylabel('y-ass')
%% 8.uzdevums
clc, clearvars, format compact, close all
y1 = @(x) cos(x.^3+3)
y2 = @(x) sin(x.^3)+3
x = -pi:0.01:pi;
plot(x,y1(x),'r-',x,y2(x),':g','LineWidth',3)
title('Trigonometriskās funkcijas')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('y=cos(x^3+3)','y2=sinx^3+3')
%% 9.uzdevums
clc, clearvars, format compact, close all
y1 = @(x) sin(x.^2+1)
y2 = @(x) cos(2*x).^2
x = -2*pi:0.01:2*pi;
zim = plot(x,y1(x),'b--',x,y2(x),'r-','LineWidth',3)
set(zim(1),'LineWidth',2)
set(zim(2),'LineWidth',3)
title('Divu funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('y=sin(x^2+1)','y=cos^2(2x)')
%% 10.uzdevums
clc, clearvars, format compact, close all
syms x
f1(x) = sqrt((cos(5*x + 3*x.^3)).^5 + (5*x.^6).^(1/7));
f2(x) = log(5*x.^3 + 5);
f3(x) = f1(x) + f2(x);

x_ver = 0.2:0.2:3.2;
y_ver = double(f1(x_ver));
tabula = [x_ver' y_ver'];
disp('x_mezgli y_mezgli')
disp(tabula)

hold on
plot(x_ver, y_ver, 'Ok', 'LineWidth', 2)
fplot(f2, [0.2, 3.2], 'r--', 'LineWidth', 2)
fplot(f3, [0.2, 3.2], 'b-', 'LineWidth', 2)
legend('y = f_1(x)', 'y = f_2(x)', 'y = f_3(x)')
grid on
hold off
%% 11.uzdevums
clc, clearvars, format compact, close all
syms x t
xt(t) = (2*(t.^2-1))/(1+t.^2)
yt(t) = (2*t*(t.^2-1))/(1+t.^3)
f(x) = sin(x)*cos(2*x)
fplot(xt,yt,[0,2*pi],'xr--','LineWidth', 1.5), hold on
fplot(f,[0,2*pi], '*g-','LineWidth', 2), hold off
title('Divu funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('x1(t)y1(t)','y=sinxcos2x')
%% 12.uzdevums
clc, clearvars, format compact, close all
y = @(x) (acot(2*x+5)).^3
f = @(x) 2*acot(x.^3)+5
x = -9:0.01:9
plot(x,y(x),':r',x,f(x),'g-','LineWidth', 3)
title('Divu funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('y(x)','f(x)')
%% 13.uzdevums
clc, clearvars, format compact, close all
y = @(x) cos(x+1).*log(x.^2-1)
f = @(x) (sin(x)+1).*(log(x.^2)-1)
x = 3:0.01:11
plot(x,y(x),'-k',x,f(x),'b--','LineWidth', 5)
title('Divu funkciju grafiki')
grid on,xlabel('x-ass'), ylabel('y-ass')
legend('y(x)','f(x)')