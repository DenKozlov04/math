%% 1.piemērs. fplot3 . Konstruēt 3D līniju

clc, clearvars, format compact, close all
% b)variants - funkcija definēta ( symbolic function )
syms t
x(t) = exp(-t/10)*sin(5*t);
y(t) = exp(-t/10)*cos(5*t);
z(t) = t;
fplot3(x,y,z,[-10,10],'-g','LineWidth',2)
xlabel('x-ass'), ylabel('y-ass'), zlabel('z-ass')
view([-3 -1 2])
% symbolic expression
figure,syms t
fplot3(exp(-t/10)*sin(5*t),exp(-t/10)*cos(5*t),t,[-10,10],'-g','LineWidth',2)
xlabel('x-ass'), ylabel('y-ass'), zlabel('z-ass')
view([-3,-1,2])