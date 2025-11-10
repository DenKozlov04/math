%% kd 2
%% 1 uzd
%%  Sastādīt funkcijas f(x) = (3x ^ 2 - x + 1)/(2x ^ 2 + 2x + 1) 
% interpolēt un aproksimēt
% vienādi attālinātu vērtību tabulu intervälä [2;10] ar soli x = 2 
% Interpolēt šo tabulu ar kubisko splainu (apzīmēsim šo funkciju ar h(x) ) 
% Aproksimet šo tabulu ar otras kartas polinomu (apzīmēsim šo funkciju ar g(x) ). 
% Atrast |h( x{1}|- g(x{1})| punktā x_{1} = 4.6
clc, clearvars, format compact, close all
% Definē funkciju f(x)
f = @(x) (3*x.^2 - x + 1) ./ (2*x.^2 + 2*x +3);

% Izveido vienādi attālinātu vērtību tabulu
intervals1 = 2;
intervals2 = 10;
solis = 2;
polinoms = 2; %kuras kārtas polinoms
x1 = 4.6; % punkts x1

x_values = intervals1:solis:intervals2; 
f_values = f(x_values);

% Aproksimē ar otras kartas polinomu g(x)
p_coefficients = polyfit(x_values, f_values, polinoms);
g = @(x) polyval(p_coefficients, x);

% Interpolē ar kubisko splainu h(x)
spline_coefficients = spline(x_values, f_values);
h = @(x) ppval(spline_coefficients, x);

% Atrast |h(x1) - g(x1)| punktā x1
difference = abs(h(x1) - g(x1))
%% 2 uzd
clc, clearvars,format compact, close all

% Definē funkciju f(x)
f = @(x) (x.^2 + x + 1) ./ (3*x.^2 + 2*x + 5);

% Izveido tabulu intervalā [0; 1.8] ar soli 0.2
x_values = 0:0.2:1.8;
f_values = f(x_values);

% Aproksimē ar formulu y(x) = a + b * sin(x) + c * cos(x) + d * sin(2*x)
ft = fittype({'1', 'sin(x)', 'cos(x)', 'sin(2*x)'});
[p, reg] = fit(x_values', f_values', ft);

% Aprēķināt starpību pēc moduļa punktā x_{0} = 1.15
x0 = 1.15;
difference = abs(f(x0) - p(x0)) 
%% 3rd sastādīt funkcijas f(x) vienādi attālinātu vērtību tabulu intervālā [1;9] ar solie deltax=1. Interpolēt ar kubisko splainu punktā x0=7.21
clc, clearvars, format compact, close all
f=@(x)(cos(x).^3).*sqrt(4.*nthroot(x.^2,3)+x.*sin(x));
xnodes=(1:1:9); ynodes=f(xnodes);
x0=7.21;
f_res=interp1(xnodes,ynodes,x0,'spline')
%% 4.uzd interpolēt ar ņūtona interpolācijas polinomu atrast interpolācijas polinoma vērtību punktā x0=4.6
clear all, clc ,close all, format compact
f=@(x)(nthroot(2+log(2+sin(x)),3)).*sin((x.^2)+3);
xnodes=[1:8];ynodes=f(xnodes);
coef=ynodes; % koeficientu 𝒂𝟏, 𝒂𝟐, … noteikšana
m = length(xnodes);
for k = 2:m
coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes (k:m) - xnodes (1:m+1-k));
end
% turpinājums
syms x
pol= coef(m); % polinoma konstuēšana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
% Ctrl+Enter
% turpinājums
polyn(x)=collect(pol); % vai → polyn(x)= expand(pol)
coefpol=sym2poly(polyn); % polinoma koeficienti (formātā double)
% Ctrl+Enter
% turpinājums
int_res=double(polyn(4.6))

%% 5 uzd aproksimēt ar 2. kārtas polinomu
clc, clearvars, format compact, close all
syms t, y = @(t)((1+cos(t))./(5+sqrt(1+sqrt(t)+t.^2))); %funkcija
inter1 = 0.4; %intervāls 1
inter2 = 1.8; %intervāls 2
solis = 0.2; %ar soli
polinomaPak = 2; %kuras pakāpes polinoms - poly4 ja 4
punkts = 1.13; %kļūda vai vērtība punktā
i = 0;
for x = inter1:solis:inter2
    i=i+1;
    f(i)=integral(y,0,x);
end
xpoints = (inter1:solis:inter2)'; ypoints = f';
p=polyfit(xpoints, ypoints,polinomaPak) %aproksimējošā polinoma koeficients ir pirmā vērtība

pol_ver=polyval(p,punkts) %polinoma vērtība punktā
f_ver=integral(y,0,punkts);
aproks_kluda=abs(pol_ver-f_ver) % Aprēķināt aproksimācijas kļūdu punktā x0 = punkts