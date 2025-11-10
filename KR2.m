% Konstruēt funkcijas f(x) = sin(3+sin(x)) * (3 + 2*arctg(3 + x^2))^(1/3) vienādi attālinātu
%vērtību tabulu intervālā [3,11] ar soli \(\Delta x = 2\).  
%Interpolēt tabulu ar **Nūtona interpolācijas polinomu** (apzīmēsim interpolācijas polinomu ar \( p(x) \)).  
%Interpolēt šo tabulu ar **kubisko splainu** (apzīmēsim splainu ar \( spl(x) \)).  
%Atrast vērtību \( |p(x_1) - spl(x_1)| \) punktā \( x_1 = 4.11 \).  
%Lai pareizi noapaļotu rezultātu, izmantot komandu \( \textit{fprintf} \).  
%Atbildi dod ar ēetriem cipariem aiz komata.

%1uzdevums
clc; clearvars; format compact;

f = @(x)sin(3+sin(x)).*(3+2*atan(3+x.^2)).^(1/3);
xnodes = 3:2:11;
ynodes = f(xnodes);
coef = ynodes;
m = length(xnodes);

for k = 2:m
    coef(k:m)=(coef(k:m)-coef(k-1:m-1))./(xnodes(k:m)-xnodes(1:m-k+1));
end

syms x
pol=coef(m);
for k=m-1:-1:1
    pol=pol*(x-xnodes(k))+coef(k);
end
p_x=matlabFunction(pol);
x1=4.11;
spl_val=interp1(xnodes,ynodes,x1,'spline');
p_val=p_x(x1);
result=abs(p_val-spl_val);
fprintf('%.4f\n',result);

%Sastādīt funkcijas f(x) = sin²(x) * √( ln(3 + cos(x)) + 5·x^(1/4)  vienādi attālinātu  
%vērtību tabulu intervālā [ 5 ; 13 ] ar soli \(\Delta x = 2\). Interpolēt šo tabulu ar **kubisko splainu**.  
%Kāda ir interpolācijas kļūdas \( |f(x_0) - f(x_0)| \) vērtība punktā \( x_0 = 10.22 \),  
%kur \( f(x_0) \) ir interpolācijas rezultāts Atbildi dod ar četriem cipariem aiz komata.

%2uzdevums
clc; clearvars; format compact;

f = @(x)(sin(x)).^2 .*sqrt(log(3+cos(x))+5*x.^(1/4));
xnodes=5:2:13;
ynodes=f(xnodes);
x0=10.22;
f_interp=interp1(xnodes,ynodes,x0,'spline');  
f_exact=f(x0);                                   
error_val=abs(f_interp-f_exact);
fprintf('%.4f\n',error_val);

%Sastādīt funkcijas f(x) = (2 + ln(1 + 2·x³))^(1/4) * (cos(e^x + 3) + 3·x) vienādi attālinātu  
%vērtību tabulu intervālā [3 ; 7] ar soli \(\Delta x = 1\)  
%Interpoļēt šo tabulu ar **Nūtona interpolācijas polinomu**.  
%Atrast interpolācijas polinoma vērtību punktā \( x_0 = 4.23 \).  
%Atbildi dod ar ēetriem cipariem aiz komata.

%3uzdevums
clc; clearvars; format compact;

f = @(x)(2+log(1+2*x.^3)).^(1/4) .*(cos(exp(x)+3)+3*x);
xnodes=3:1:7;
ynodes=f(xnodes);
coef = ynodes;
m=length(xnodes);
for k=2:m
    coef(k:m)=(coef(k:m)-coef(k-1:m-1))./(xnodes(k:m)-xnodes(1:m-k+1));
end

syms x
pol=coef(m);
for k=m-1:-1:1
    pol=pol*(x-xnodes(k))+coef(k);
end
p_x = matlabFunction(pol);

x0=4.23;
result=p_x(x0);
fprintf('%.4f\n',result);

%Sastādīt funkcijas f(x) = (2 + ln(1 + 2·x³))^(1/4) vienādi attālinātu  
%vērtību tabulu intervālā [2 ; 8] ar soli \(\Delta x = 1\).  
%Atrast pirmās kārtas izdalīto starpību \( f[x_1, x_2] \).  
%Atbildi dod ar ēetriem cipariem aiz komata.

%4uzdevums
clc; clearvars; format compact;
f = @(x)(2+log(1+2*x.^3)).^(1/4);
xnodes=2:1:8;
ynodes=f(xnodes);

f_x1_x2=(ynodes(2)-ynodes(1))/(xnodes(2)-xnodes(1));
fprintf('%.4f\n', f_x1_x2);

%Sastādīt funkcijas f(x) = (x⁴)^(1/5) + (2·x + x³) * sin(x + 6) vienādi attālinātu  
%vērtību tabulu intervālā [1 ; 13] ar soli \(\Delta x = 2\).  
%Interpolēt šo tabulu ar Nūtona interpolācijas polinomu.  
%Uzzīmēt Nūtona interpolācijas polinoma grafiku dotajā intervālā.

%5uzdevum
clc; clearvars; format compact;

f = @(x)(x.^4).^(1/5)+(2*x+x.^3).* sin(x+6);
xnodes=1:2:13;
ynodes=f(xnodes);
coef=ynodes;
m=length(xnodes);
for k=2:m
    coef(k:m)=(coef(k:m)-coef(k-1:m-1))./(xnodes(k:m)-xnodes(1:m-k+1));
end

syms x
pol=coef(m);
for k=m-1:-1:1
    pol=pol*(x-xnodes(k))+coef(k);
end
p_x=matlabFunction(pol);

x_plot=1:0.1:13;
y_pol=p_x(x_plot);
y_exact=f(x_plot);

figure;
plot(x_plot,y_pol,'r-','LineWidth',2);
hold on;
plot(xnodes,ynodes,'bo','MarkerSize',8,'MarkerFaceColor','b');
plot(x_plot, y_exact,'g--','LineWidth',1);
grid on;

%Вторая часть 4 домашки

%Sastādīt funkcijas f(x) = x^(4/5) + (2·x + x²) * ln(x) vienādi attālinātu  
%vērtību tabulu intervālā [1 ; 13] ar soli \(\Delta x = 2\).  
%Interpolēt šo tabulu ar **Nūtona interpolācijas polinomu** (\( pol(x) \)).  
%Uzzīmēt interpolācijas kļūdas \( |f(x) - pol(x)| \) grafiku.

%1uzdevums
clc;clearvars;format compact;

f=@(x)x.^(4/5)+(2.*x+x.^2).*log(x);
xnodes=1:2:13;
ynodes=f(xnodes);
coef=ynodes;
m=length(xnodes);

for k=2:m
    coef(k:m)=(coef(k:m)-coef(k-1:m-1))./(xnodes(k:m)-xnodes(1:m-k+1));
end

syms x
pol=coef(m);
for k=m-1:-1:1
    pol=pol*(x-xnodes(k))+coef(k);
end
p_x=matlabFunction(pol);

xq=linspace(1,13,500);
err=abs(f(xq)-p_x(xq));

figure;
plot(xq,err,'LineWidth',1.2);grid on;

%Konstruēt funkcijas f(x) = ln(5 + cos(2·x) + sin(x)) * (5 + arctan(x))^(1/4) vienādi attālinātu  
%vērtību tabulu intervālā [2,6] ar soli \(\Delta x = 1\).  
%Interpolēt šo tabulu ar kubisko splainu. Uzzīmēt splaina grafiku dotajā intervālā.

%2uzdevums
clc;clearvars;format compact;

f=@(x)log(5+cos(2.*x)+sin(x)).*(5+atan(x)).^(1/4);
xnodes=2:1:6;
ynodes=f(xnodes);
xq=linspace(2,6,500);
spl=interp1(xnodes,ynodes,xq,'spline');

figure;
plot(xq,spl,'LineWidth',1.2);grid on;

%Sastādīt funkcijas f(x) = (2 + ln(1 + 2·x³))^(1/4) vienādi attālinātu  
%vērtību tabulu intervālā [2 ; 8] ar soli \(\Delta x = 1\).  
%Interpolēt šo tabulu ar **Nūtona interpolācijas polinomu**.  
%Atrast interpolācijas polinoma koeficientu pie \( x^4 \).  
%Atbildi dod ar četriem cipariem aiz komata.

%3uzdevums
clc;clearvars;format compact;

f=@(x)(2+log(1+2.*x.^3)).^(1/4);
xnodes=2:1:8;
ynodes=f(xnodes);
coef=ynodes;
m=length(xnodes);

for k=2:m
    coef(k:m)=(coef(k:m)-coef(k-1:m-1))./(xnodes(k:m)-xnodes(1:m-k+1));
end

syms x
pol=coef(m);
for k=m-1:-1:1
    pol=pol*(x-xnodes(k))+coef(k);
end

p=poly2sym(sym2poly(pol),x);
a4=double(coeffs(p,x,'All'));
deg=length(a4)-1;
idx=deg-4+1;
fprintf('%.4f\n',a4(idx));

%Konstruēt funkcijas f(x) = arccot(3 + 3·x + x²) * sin(3 + ln(3 + x)) vienādi attālinātu  
%vērtību tabulu intervālā [2,10] ar soli \(\Delta x = 2\).  
%Interpolēt šo tabulu ar kubisko splainu (\( spl(x) \)).  
%Uzzīmēt interpolācijas kļūdas \( | f(x) - spl(x) | \) grafiku.

%4uzdevums
clc;clearvars;format compact;

f=@(x)acot(3+3.*x+x.^2).*sin(3+log(3+x));
xnodes=2:2:10;
ynodes=f(xnodes);
xq=linspace(2,10,1000);
spl=interp1(xnodes,ynodes,xq,'spline');
err=abs(f(xq)-spl);

figure;
plot(xq,err,'LineWidth',1.2);grid on;

%Sastādīt funkcijas f(x) = cos(2 + sin(3·x)) * (3 + 2·ln(1 + 3·x))^(1/3) vienādi attālinātu  
%vērtību tabulu intervālā [4 : 16] ar soli \(\Delta x = 3\).  
%Interpolēt šo tabulu ar **kubisko splainu**. Atrast splaina vērtību punktā \( x_0 = 8.43 \).  
%Atbildi dod ar četriem cipariem aiz komata.

%5uzdevums
clc;clearvars;format compact;

f=@(x)cos(2+sin(3.*x)).*(3+2.*log(1+3.*x)).^(1/3);
xnodes=4:3:16;
ynodes=f(xnodes);
x0=8.43;
spl_val=interp1(xnodes,ynodes,x0,'spline');
fprintf('%.4f\n',spl_val);


%5 дз 1 часть

%Sastādīt funkcijas f(x_i) = ∫₀^x_i cos(12·t² + √(t³ + 6)) dt vienādi attālinātu vērtību tabulu intervālā [1.4 : 2] ar soli 0.1.
%Interpolēt tabulu ar Nūtona interpolācijas polinomu \( p_n(x) \).
%Aprēķināt interpolācijas kļūdu \( |f(x_0) - p_n(x_0)| \) punktā \( x_0 = 1.74 \).
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums1
clear, clc
a =1.4; 
b =2; 
h =0.1;
x =a:h:b;
n =length(x);
for i =1:n
    f(i) =integral(@(t) cos(12*t.^2+sqrt(t.^3 + 6)),0,x(i));
end

d =zeros(n,n);
d(:,1) =f';
for j =2:n
    for i =1:(n-j+1)
        d(i,j) =(d(i+1,j-1)-d(i,j-1))/(x(i+j-1)-x(i));
    end
end

x0 =1.74;
p =d(1,1);
q =1;
for j =2:n
    q =q*(x0-x(j-1));
    p =p+d(1,j)*q;
end

f0 =integral(@(t)cos(12*t.^2+sqrt(t.^3+6)),0,x0);
k =abs(f0-p);
fprintf('%.4f\n',x0,x0,k);


%Sastādīt funkcijas f(x_i) = ∛(4 + x_i² + x_i³) ÷ (3 + cos(x_i) + √(2 + sin(2·x_i))) vienādi attālinātu vērtību tabulu intervālā [1:9] ar soli 1.
%Aproksimēt tabulu ar formulu \( y(x) = a + 2bx + c \cos 2x \), izmantojot komandu fit. Uzzīmēt aproksimējošās funkcijas grafiku dotajā intervālā.

%uzdevums2
clear, clc
a =1; b =9; h =1;
x = a:h:b;
n =length(x);
for i =1:n
    f(i) =(4+x(i)^2+x(i)^3)^(1/3)/(3+cos(x(i))+sqrt(2+sin(2*x(i))));
end

for i = 1:n
    A(i,1) =1;
    A(i,2) =2*x(i);
    A(i,3) =cos(2*x(i));
end

coef =A\f';
a1 =coef(1); 
b1 =coef(2); 
c1 =coef(3);
xx =a:0.01:b;
for i =1:length(xx)
    y(i) =a1+2*b1*xx(i)+c1*cos(2*xx(i));
end

plot(x,f,'bo','MarkerFaceColor','b');hold on
plot(xx,y,'r-','LineWidth',1.5)
grid on

%Sastādīt funkcijas f(x_i) = ∫₀^x_i (1 + √(t³)) ÷ (5 + cos(t) + √(5 + sin(2·t) + t²)) dt vienādi attālinātu vērtību tabulu intervālā [5: 14] ar soli 1.
%Aproksimēt šo tabulu ar 4 pakāpes polinomu, izmantojot komandu polyfit.  
%Kāds ir aproksimējošā polinoma koeficients pie \( x^3 \)?
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums3
clear, clc
a = 5;
b = 14;
h = 1;
x =a:h:b;
n =length(x);
for i =1:n
    f(i) =integral(@(t)(1 + sqrt(t.^3))./(5+cos(t)+sqrt(5+sin(2*t)+t.^2)),0,x(i));
end
p = polyfit(x,f,4);
fprintf('%.4f\n',p(2));

%Sastādīt funkcijas f(x_i) = (x_i² + cos(x_i)) ÷ (1 + 3·x_i + √(x_i² + 1)) vienādi attālinātu vērtību tabulu intervālā \[ [0.6 : 2.4] \] ar soli 0.2.
%Aproksimēt šo tabulu ar 4 pakāpes polinomu, izmantojot komandu polyfit.
%Kāds ir aproksimējošā polinoma koeficients pie \( x^3 \)?
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums4
clear, clc
a =0.6;
b =2.4;
h =0.2;
x =a:h:b;
n =length(x);
for i = 1:n
    f(i) = (x(i)^2+cos(x(i)))/(1+3*x(i)+sqrt(x(i)^2+1));
end
p =polyfit(x,f,4);
fprintf('%.4f\n',p(2));

%Sastādīt funkcijas f(x_i) = (x_i² + ln(2·x_i + 5)) ÷ (x_i² + √(2·x_i + 3) + 2) vienādi attālinātu vērtību tabulu intervālā
%[1:10] ar soli \(\Delta x = 1\)
%Aproksimēt šo tabulu ar 5 pakāpes polinomu, izmantojot komandu polyfit.
%Kāda ir šī polinoma vērtība punktā \( x_0 = 5.339 \)
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums5
clear, clc
a =1;
b =10;
h =1;
x =a:h:b;
n =length(x);
for i =1:n
    f(i) =(x(i)^2+log(2*x(i)+5))/(x(i)^2+sqrt(2*x(i)+3)+2);
end

p =polyfit(x,f,5);
x0 =5.33;
y0 =polyval(p,x0);
fprintf('%.4f\n',x0,y0);


%5 дз 2 часть

%Sastādīt funkcijas f(x_i) = ∛(4 + x_i² + x_i³) ÷ (3 + cos(x_i) + √(2 + sin(2·x_i))) vienādi attālinātu vērtību tabulu intervālā [1 : 9] ar soli 1.
%Aproksimēt tabulu ar formulu \( y(x) = a + 2bx + c\cos 2x \), izmantojot komandu fit.
%Kāda ir aproksimācijas kļūdas \( |y(x_0) - f(x_0)| \) vērtība punktā \( x_0 = 5.67 \)
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums1
clear, clc
a =1;
b =9;
h =1;
x =a:h:b;
n =length(x);

for i =1:n
    f(i) =(4+x(i)^2+x(i)^3)^(1/3)/(3+cos(x(i))+sqrt(2+sin(2*x(i))));
end
for i =1:n
    A(i,1) =1;
    A(i,2) =2*x(i);
    A(i,3) =cos(2*x(i));
end
coef =A\f';
a1 =coef(1); 
b1 =coef(2);
c1 =coef(3);

x0 =5.67;
y0 =a1+2*b1*x0+c1*cos(2*x0);
f0 =(4+x0^2+x0^3)^(1/3)/(3+cos(x0)+sqrt(2+sin(2*x0)));
k =abs(y0-f0);
fprintf('%.4f\n',k);


%Sastādīt funkcijas f(x_i) = 1 ÷ ((1 + x_i)² + √(ln(x_i + 2) + 3·x_i)) vienādi attālinātu vērtību tabulu intervālā [1.2; 3] ar soli \(\Delta x = 0.2\).
%Aproksimēt šo tabulu ar formulu \( y(x) = a + bx + c \cos x \) izmantojot komandu fit.
%Aprēķināt \( y(x_0) \), ja \( x_0 = 1.33 \). Atbildi dod ar četriem cipariem aiz komata.

%uzdevums2
clear, clc
a =1.2;
b =3;
h =0.2;
x =a:h:b; 
n =length(x);

for i =1:n
    f(i) =1/((1+x(i))^2+sqrt(log(x(i)+2)+3*x(i)));
end

for i =1:n
    A(i,1) =1;
    A(i,2) =x(i);
    A(i,3) =cos(x(i));
end
coef =A\f';
a1 =coef(1); 
b1 =coef(2); 
c1 =coef(3);

x0 =1.33;
y0 =a1+b1*x0+c1*cos(x0);
fprintf('%.4f\n',x0,y0);

%Sastādīt funkcijas f(x_i) = ∫₀^x_i (1 + √(t³)) ÷ (5 + cos(t) + √(5 + sin(2·t) + t²)) dt vienādi attālinātu vērtību tabulu intervālā [5: 14] ar soli 1.
%Aproksimēt šo tabulu ar 4.pakāpes polinomu, izmantojot komandu polyfit.  
%Kāda ir šī polinoma vērtība punktā \( x_0 = 8.88 \)  
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums3
clear, clc
a =5;
b =14;
h =1;
x =a:h:b;
n =length(x);

for i =1:n
    f(i) =integral(@(t)(1+sqrt(t.^3))./(5+cos(t)+sqrt(5+sin(2*t)+t.^2)),0,x(i));
end
p =polyfit(x,f,4);
x0 =8.88;
y0 =polyval(p,x0);
fprintf('%.4f\n',x0,y0);

%Sastādīt funkcijas f(x_i) = ∫₀^x_i (exp(√t · sin(3·t² + 1)) + cos(2·t)) dt vienādi attālinātu  
%vērtību tabulu intervālā [3: 9] ar soli 1.  
%Interpolēt tabulu ar kubisko splainu \( s(x) \). Aprēķināt \( s(3.42) \).  
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums4
clear, clc
a =3;
b =9;
h =1;
x =a:h:b;
n =length(x);

for i =1:n
    f(i) =integral(@(t)exp(sqrt(t).*sin(3*t.^2+1))+cos(2*t),0,x(i));
end
s =spline(x,f);
x0 =3.42;
y0 =ppval(s,x0);
fprintf('%.4f\n',x0,y0);


%Sastādīt funkcijas f(x_i) = ∛(4 + x_i² + x_i³) ÷ (3 + cos(x_i) + √(2 + sin(2·x_i))) vienādi attālinātu vērtību tabulu intervālā [1:9] ar soli 1.
%Interpolēt šo tabulu ar kubisko splainu (apzīmēsim šo funkciju ar \( h(x) \)).  
%Aproksimēt šo tabulu ar 5.kārtas polinomu, izmantojot komandu polyfit (apzīmēsim šo funkciju ar \( g(x) \)).  
%Atrast \( |h(x_0) - g(x_0)| \) punktā \( x_0 = 2.44 \).  
%Atbildi dod ar četriem cipariem aiz komata.

%uzdevums5
clear, clc
a =1;
b =9;
h =1;
x =a:h:b;
n =length(x);

for i =1:n
    f(i) = (4+x(i)^2+x(i)^3)^(1/3)/(3 +cos(x(i))+sqrt(2+sin(2*x(i))));
end
hs =spline(x,f);
p =polyfit(x,f,5);

x0 =2.44;
h0 =ppval(hs,x0);
g0 =polyval(p,x0);
k =abs(h0-g0);
fprintf('%.4f\n',k);


