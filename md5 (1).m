%% 1. uzd % Aproksimē tabulu pēc prasītās formulas ar fit, tad atrod grafiku kuru vajag.
% Šajā uzdevuma prasītā formula:
% y(x) = a + b*sqrt(x) + c*(sin(x) + 5)
clc, clearvars, format compact, close all
f = @(x) ( x .* exp(-0.1.*x) )./( 2 + cos(2.*x) + sqrt( 3 + sin(2.*x) ) )
xpoints = (5:1:13)'; ypoints = f(xpoints);

ft = fittype({'1','sqrt(x)','sin(x) + 5'});
[y,reg] = fit(xpoints,ypoints,ft)
figure
plot(y,xpoints,ypoints), grid on
title('Aproksimācija ar trig. funkciju y = a+b*sqrt(x)+c*(sinx + 5)')

%% 2. uzd % Sastādīt tabulu integrālim no f(x) jeb kodā y(x), tad aproksimēt ar 4. kārtas polinomu izmantojot polyfit
% Un atrast koeficientu pie x0 vai pēc uzd. prasībām (laikam)
clc, clearvars, format compact, close all
syms t, y = @(t) ( 1 + sin(t.^3) )./( 5 + cos(4.*t) + 2.*t.^2 );

i = 0;
for x = 0.4:0.2:2.2 % tabula integrālim no f(x) intervālā [0.4;2.2] ar soli 0.2
    i = i+1;
    f(i) = integral(y,0,x);
end

xpoints = (0.4:0.2:2.2)'; ypoints = f';
m = length(xpoints);
p = polyfit(xpoints,ypoints,3);
xapprox = xpoints(1):0.01:xpoints(m); yapprox = polyval(p,xapprox);
fun_prob1(p)

%% 3. uzd % Aproksimēt funkciju f(x) ar 4. kārtas polinomu izmantojot polyfit
% Atrast koeficientu pie x0 vai pēc uzd. prasībām
clc, clearvars, format compact, close all
f = @(x)(x.^2 + log(2.*x+5))./(x.^2 + sqrt(2.*x+3) + 2); %nthroot ir n-tā sakne
xpoints = (3:1:12)'; ypoints = f(xpoints); % [3,12] solis=1
p = polyfit(xpoints,ypoints,4); % Aproksimācija pēc polyfit
fun_prob1(p) % Te atrod koeficientus, KODA LĒJĀ !!!-->ĀREJĀ FUNKCIJA<--!!!

%% 4. uzd % Aproksimācija ar 4. kārtas polinomu (Atrast p vērtību pie x0 = 1.33)
clc, clearvars, format compact, close all
f = @(x)1./( (1+x).^2 + sqrt(log(x+2) + 3.*x) );
xpoints = (1.2:0.2:3)'; 
ypoints = f(xpoints);

p = polyfit(xpoints,ypoints,4); % aprokismācija ar 4. kārtas pol
x0 = 1.33;
polinoma_vertiba_x0 = polyval(p,x0)

%% 5. uzd %% Sastādīt tabulu integrālim no f(x) jeb kodā y(x), tad interpolēt ar Ņūtona polinōmu
% Pēc tam atrast interpolācijas kļūdu abs(f_x0 - newton_x0)
clc, clearvars, format compact, close all
syms t, y = @(t)(exp(sqrt(t) .* sin(3.*t.^2 + 1)) + cos(2.*t));

i = 0;
for x = 3:1:9 % tabula f(x) intervālā [3;9] ar soli 1
    i = i+1;
    f(i) = integral(y,0,x);
end

xnodes = (3:1:9)'; ynodes = f'; % tabula f(x) intervālā [3;9] ar soli 1
m = length(xnodes)

coef = ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./...
               (xnodes(k:m) - xnodes(1:m+1-k));
end

syms x, pol = coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol = pol*(x-xnodes(k))+coef(k);
end

polyn(x) = collect(pol);
coefpol = sym2poly(polyn); % polinoma koeficienti
fun_prob6(coefpol);        % polinoma drukâðana


x0 = 4.44;
newton_x0 = double(polyn(x0)); % Ņūtona polinoma vērtība x0 = 4.44
f_x0 = integral(y,0,x0); % funkcijas vērtība punktā x0 = 4.44
interpolacijas_kluda = abs(f_x0 - newton_x0);
fprintf('Atbilde: %.4f\n', interpolacijas_kluda); % Ar 4 zīmīgajiem cipariem aiz komata

%% ārējas funkcijas %%  ārējas funkcijas  %%   ārējas funkcijas %%%%%%%%%%
% pamatprogrammas beigas

%% Otrā daļa
%% 1. uzd % Aproksimē tabulu pēc prasītās formulas ar fit
% Šajā uzdevumā prasītā formula:
% y(x) = a + b*sqrt(x) + c*sin(x)
% Pēc tam atrod vērtību pie x0 pēc šīs formulas
clc, clearvars, format compact, close all
f = @(x)( 1 + 3.*x + x.^2 ) ./ ( 1 + 2.*x + sqrt(3 + sin(x)));
xpoints = (3:1:12)'; ypoints = f(xpoints);

ft = fittype({'1','sqrt(x)','sin(x)'});
[y,reg] = fit(xpoints,ypoints,ft)

x0 = 6.33
y_pie_x0 = y(x0)

%% 2. uzd % Sastāda vērtību tabulu pēc integrāļa no f(x) vai kodā y(x)
% Tad interpolē ar kubisko splainu, atrod splaina vērt. pie x0 = 3.42
clc, clearvars, format compact, close all
syms t, y = @(t)(exp(sqrt(t) .* sin(3.*t.^2 + 1)) + cos(2.*t));

i = 0;
for x = 3:1:9 % tabula integrālim no f(x) % intervāls: [3;9] ar soli 1
    i = i+1;
    f(i) = integral(y,0,x);
end
xpoints = (3:1:9); ypoints = f % intervāls: [3;9] ar soli 1

x0 = 3.42; % s(3.42)

h_x0 = interp1(xpoints,ypoints,x0,'spline') % (h(x0))

%% 4. uzd  Sastādīt tabulu integrālim no f(x) jeb kodā y(x)
% Tad aproksimē šo tabulu ar 3. kārtas polinomu, aprēķina polinoma vērtību pie x0 = 1.76
clc, clearvars, format compact, close all
syms t, y = @(t)( 3 + 3.*t.^2 ) ./ ( sqrt(2.*t + log(4 + 3.*t + 2.*t.^2) ) );

i = 0;
for x = 0:0.5:4.5 % tabula integrālim no f(x) intervālā [0;4.5] ar soli 0.5
    i = i+1;
    f(i) = integral(y,0,x);
end

xpoints = (0:0.5:4.5)'; ypoints = f'; % intervālā [0;4.5] ar soli 0.5
m = length(xpoints);
p = polyfit(xpoints,ypoints,4)

x0 = 1.76;
vertiba_pie_x0 = polyval(p,x0)

%% 5. uzd %% Aproksimē tabulu pēc prasītās formulas ar fit
% Pēc tam atrod vērtību pie x0 pēc šīs formulas un funkcijas vērt. pie x0
clc, clearvars, format compact, close all
f = @(x) ( 3 + log(2.*x + 3) ) ./ ( 2 + sqrt(3 + sin(3 + 3.*x + x.^2)) );
xpoints = (3:0.2:4.8)'; ypoints = f(xpoints); % intervālā [3;4.8] ar soli 0.2
ft = fittype({'1','x','cos(2*x)'});
[y,reg] = fit(xpoints,ypoints,ft)

x0 = 3.73; % vērtība punktā x0 = 3.73
apr_error = abs(y(x0)-f(x0))

% ârçja funkcija (1.piemçrs). Aproksimâcija
% polinoma drukâðana 
function fun_prob1(koef) 
   m = length(koef);
   fprintf('\n Atbilde.  %.0f.kârtas polinoms: \n  ',m-1)
   n = m-1;
   for i = 1:m
      if koef(i) < 0
        fprintf(' %.4fx^%.0f',koef(i),n)
      else
         fprintf(' +')
         fprintf('%.4fx^%.0f',koef(i),n) 
      end
      n = n-1;
   end
   fprintf('\n')
end

% ârçja funkcija (6.piemçrs). Aproksimâcija
% polinoma drukâðana 
function fun_prob6(koef) 
   m = length(koef);
   fprintf('\n Atbilde.  %.0f.kârtas polinoms: \n  ',m-1)
   n = m-1;
   for i = 1:m
      if koef(i) < 0
        fprintf(' %.4fx^%.0f',koef(i),n)
      else
         fprintf(' +')
         fprintf('%.4fx^%.0f',koef(i),n) 
      end
      n = n-1;
   end
   fprintf('\n')
end