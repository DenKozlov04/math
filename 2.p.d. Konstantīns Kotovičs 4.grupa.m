%% Konstantīns Kotovičs 4.grupa
%% 1.uzd.
clc, clearvars,format compact, close all
f = @(x) (x.^3.*log(x))./(exp(x)+x+sqrt(1+sin(x)));
xpoints = (0.2:0.2:1.2)' 
ypoints = f(xpoints)
plot(xpoints,ypoints,'or','LineWidth',3), grid on
title('Datu punkti')

% turpinвjums – aproksimвcija ar otras kвrtu polinomu 
p = polyfit(xpoints,ypoints,2)
xapprox = 0.2:0.01:1.2; yapprox = polyval(p,xapprox);
figure
plot(xpoints,ypoints,'or',xapprox,yapprox,'k','LineWidth',3)
grid on
title('Aproksimвcija ar 2.kвrtas polinomu')

x0 = 1.03; % vērtība punktā x0 = 1.03

apr_error = abs(polyval(p,x0)-f(x0))



%% 2. uzd.
clc, clearvars,format compact, close all
syms t, y = @(t) (3+log(6+t.*cos(t)))./(1+sqrt(atan(2.*t+4)+2));
i = 0;
for x = 0.6:0.2:1.8
    i = i+1;
    f(i) = integral(y,0,x);
end
xpoints = (0.6:0.2:1.8)'; ypoints = f';
p = polyfit(xpoints,ypoints,3)

% polinoma vзrtоba punktв x1=1.84  (4)
x1 = 1.84; pol_value = polyval(p,x1)
fprintf('%.4f\n', pol_value)

%% 3.uzd.
clc, clearvars, format compact, close all 
f = @(x) sin(sqrt(x)+1) .* sqrt(cos(1+x) + log(x) + 3);
xnodes = (1:2:11); ynodes=f(xnodes);
x1 = 6.3;
int_res=interp1(xnodes,ynodes,x1,'spline')

% turpinājums
disp('Atbilde:')
fprintf('interpolācijas rezultāts punktā(6.3) = %.4f\n', int_res)

%% 4. uzd.
clc, clearvars, format compact, close all 
syms x, f = @(x) x.^2 .* ((x).^(4/3) + cos(x)).^(1/2);

xnodes = (1:2:9); ynodes = f(xnodes);  coef = ynodes;
m = length(xnodes);
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./...
                (xnodes(k:m) - xnodes(1:m+1-k));
end
pol= coef(m); % polinoma konstuēšana
for k = m-1:-1:1
   pol = pol*(x-xnodes(k))+coef(k);
end

polyn(x)=collect(pol)   % simboliskā funkcija
coefpol=sym2poly(polyn) % polinoma koeficienti

int_res = double(polyn(6.4))   % polinoma vērtība punktā x0=6.4

% turpinājums
disp('Atbilde:')
fprintf('Interpolācijas rezultāts punktā(6.4) = %.4f\n', int_res)


%% 5.uzd.
clc, clearvars,format compact, close all
f = @(x) atan(2+x.^2+3.*x.^3);
xpoints = (1:2:9)'
ypoints = f(xpoints)
plot(xpoints,ypoints,'or','LineWidth',3)
grid on
title('Datu punkti')

% turpinвjums – aproksimвcija ar 2.kвrtas polinomu
format longG
g = polyfit(xpoints,ypoints,2)

m = length(xpoints);
xapprox = xpoints(1):0.01:xpoints(m); yapprox = polyval(g,xapprox);
figure
plot(xpoints,ypoints,'or',xapprox,yapprox,'k','LineWidth',3)
grid on
title('Aproksimвcija ar 2.kвrtas polinomu')

% turpinвjums – interpolвcija ar kubisko splainu punktв x0
x1 = 4.3;
h_x1 = interp1(xpoints,ypoints,x1,'spline')

% turpinвjums. Aprзнinвt kпыdu punktв x_1 = 4.3
% (aproksimвcija un splaina)
%( kв starpоbu pзc moduпa starp funkciju h_x1 un g(x)
% vзrtоbвm punktв x_1). 
kluda = abs(h_x1 - polyval(g,x1))
fprintf('Atbilde: \n ')
fprintf('  Kпыda punktв (4.3) = %.4f \n ',kluda)
format