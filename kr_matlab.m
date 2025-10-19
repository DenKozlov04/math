%% 1.
% Dota sistema AX=B, kur A=(4x4). B=(4x1).
% Pielietot 'vienkarso iteraciju metodi' un izpildit '8 iteracijas',
% izmantojot sakuma tuvinajumu X^(0) = (0 0 0 0)^T.
% Izmantot iteraciju parametru tau=0.01. Kada ir tuvinata atrisinajuma norma
% ||X^(8)||2 8.iteracija? Atbildi dod ar cetriem cipariiem aiz komata.

clc; clearvars; format compact;

A = [6,4,8,10;4,12,12,4;8,12,144,9;10,4,9,36];
b = [3;0;5;7];
tau = 0.01; % tau

% x = [1;-2;0.5;3]; если x не 0
x = zeros(4,1);
k = 0; %или x = [0;0;0,0;];

while k < 8
    x = x + tau*(b - A*x);  % vienkāršā iterācija
    k = k + 1;
end

norm_x8 = norm(x,2);
fprintf('||X^(8)||_2 = %.4f\n', norm_x8);

%% 2.
% Dota sistema AX=B, kur A=(4x4), B=(4x1).
% Pienemot, ka sistēmu risina ar 'vienkarso iteraciju metodi',
% noskaidrot, par kadu skaitli ir jabut mazakai iteraciju parametra tau
% vertibu, lai metode konvergetu? Atbildi dod cetriem cipariem aiz komata.

clc; clearvars; format compact;

A = [4,2,4,6;2,7,8,2;4,8,49,5;6,2,5,21];
b = [1;7;5;4]; %#ok<NASGU>  % pievienots pilnīgumam, uz τ neietekmē

lam_max = max(eig(A));   % A ir simetriska
tau_limit = 2/lam_max;   % nepieciešams tau < tau_limit
fprintf('tau jabut < %.4f\n', tau_limit);

%% 3.
% Dota vienadojumu sistema AX=B, kur A=(3x4), X=(x1,x2,x3,x4), B=(3x1).
% Atrisinat so sistemu 'ar Gausa metodi', izmantojot komandu 'rref'.
% Atrast x3, ja x4 =7

clc; clearvars; format compact;

A = [1,3,2,4;2,-1,5,2;-1,3,2,1];
b = [-1;3;5];
x4 = 7; % x2 = 5

% Pārvietojam 4. kolonnas ieguldījumu uz labo pusi un risinām 3x3 sistēmu
A3 = A(:,1:3); % A3 = A(:, [1 3 4]);
b3 = b - A(:,4)*x4; %b3 = b - A(:,2)*x2;

R = rref([A3 b3]);     % Gausa metode (reduced row-echelon form)
x = R(:,end);          % [x1; x2; x3]
x3 = x(3); %x1 = x(1); x3 = x(2); x4 = x(3);
fprintf('x3 = %.4f\n', x3);

%% 4. 
% Dota vienadojumu sistema AX=B, kur A=(3x3), B(3x1).
% Sakuma tuvinajums X^(0) =(000)^T.
% Ar 'Jakobi metodi' izpildit '9 iteracijas'.
% Cik liela ir norma ||e^(7)||2 7.iteracija, kur e^(n) = X^(n) - X^(n-1)?

clc; clearvars; format compact;

A = [1,3,4;2,5,1;2,-3,7];
b = [1;1;9];

% Jacobi: x^{k+1} = D^{-1}(b - (L+U)x^{k})
D = diag(diag(A));
R = A - D;
Dinv = inv(D);

x = zeros(3,1);
eps_list = cell(9,1);

for k = 1:9
    x_new = Dinv*(b - R*x);
    eps_list{k} = x_new - x;   % ε^(k)
    x = x_new;
end

eps7 = eps_list{7};
val = norm(eps7,2);
fprintf('||ε^(7)||_2 = %.4f\n', val);

%% 5.
% Dota vienadojumu sistema AX=B, kur A=(4x3), X=(4x1), B=(3x1).
% Atrisinat so sistemu 'ar Gausa metodi', izmantojot komandu 'rref'.
% Atrast visparigo atrisinajumu.
clc; clearvars; format compact;
syms x1 x2 x3 x4

A = [4,3,5,1;2,1,-3,2;1,4,4,1];
b = [1;7;9];

eqs = A*[x1;x2;x3;x4] == b;
S = solve(eqs, [x1 x2 x3]);

fprintf('x1 = %s\n', char(S.x1));
fprintf('x2 = %s\n', char(S.x2));
fprintf('x3 = %s\n', char(S.x3));


%% 6.
% Dota vienadojumu sistema AX=B, kur A=(3x4), B=(3x1).
% Sakuma tuvinajums X^(0) =(000)^T.
% Ar 'Jakobi metodi' izpildit '8 iteracijas'
% Kada ir atrisinajuma norma||X^(8)||2 8.iteracija?
% Atbildi dod ar cetriem cipariem aiz komata
clc; clearvars; format compact;

A = [6,1,3;2,-4,7;3,-1,5];
b = [-3;2;0];

% Jacobi: x^{k+1} = D^{-1}(b - (L+U)x^{k})
D = diag(diag(A));
R = A-D;
x = zeros(3,1);

k = 0;
while k<8
    x = D\(b-R*x);  % Jakobi iterācija
    k = k+1;
end

val = norm(x,2);
fprintf('||X^(8)||_2 = %.4f\n',val);
%=================================================================
%мб будет
%
%Dota vienādojumu sistēma AX-B, kur A =...
%Izpildīt 15 iterācijas ar vienkāršo iterāciju metodi, 
%izmantojot sākuma tuvinājumu X(0) = (0 0 0 0)T".
%Aprēķinos izmantot optimālo parametra tau vērtību. 
%Atrisināt šo sistēmu arī ar komandu linsolve 
%(ar komandu linsolve iegūto rezultātu apzīmēsim ar Xs). 
%Cik liela ir norma ||Xs - X(15)||? 
%Atbildi dod ar četriem cipariem aiz komata.

A = [7,1,3,8;1,9,14,2;3,14,81,4;8,2,4,27];
B = [-7.5;2.3;-6.7;8.9];
x = zeros(4,1);

% 2. opt solis
lambda = eig(A);
tau_opt = 2/(max(lambda)+min(lambda));

% 15 iteracijas
for k = 1:15
    x = x - tau_opt*(A*x - B);
end
x15 = x;

% 4. Точное решение с использованием linsolve
x_s = linsolve(A,B);

norma = norm(x_s - x15);
fprintf('%.4f\n', norma);

%Dota vienādojumu sistēma AX=B, kur A=...
%Pieņemot, ka sistēmu risina ar vienkāršo iterāciju metodi, 
% noskaidrot, par kādu skaitli ir jābūt mazākai iterāciju parametra tau vērtībai, 
%lai metode konverģētu?
%Atbildi dod ar četriem cipariem aiz komata.
A = [28,2.5,6.5,10.1;2.5,13,16,4;6.5,16,144,6;10.1,4,6,36];
B = [3;-12;3;-4.8];
x = zeros(4,1);

lambda = eig(A);
tau_max = 2/max(lambda);

fprintf('konveerge, ja 0 < tau < %.4f\n', tau_max);

%
%Dota vienādojumu sistēma AX B, kur A
%Pielietot minimālās nesaistes metodi x^(n+1) = x^(n) + Taun+1(B − Ax^(n)),
%n = 0,1,... 
%un izpildīt 19 iterācijas, izmantojot sākuma tuvinājumu X(0) = (-1 -1 -1)T
%Atrisināt šo sistēmu arī ar komandu linsolve
%(ar komandu linsolve iegūto rezultātu apzīmēsim ar Xs).
%Cik liela ir norma ||Xs - X(19)||2? Atbildi dod ar četriem cipariem aiz komata.

A = [5.3,6.4,8;6.4,25.3,4.6;8,4.6,31];
B = [9.6;,7.6;,3.6];
x = [-1;-1;-1];

for k = 1:19
    r = B-A*x;
    Ar = A*r;
    tau = (r'*Ar)/(Ar'*Ar);
    x = x+tau*r;
end
x19 = x;
x_s = A\B;
norma = norm(x_s-x19);
fprintf('%.4f\n',norma);

%Izpildīt 29 iterācijas ar vienkāršo iterāciju metodi, 
% izmantojot sākuma tuvinājumu X(0) = (0 0 0 0)^T. 
% Aprēķinos izmantot optimālo parametra Tau vērtību. 
% Atrisināt šo sistēmu arī ar komandu linsolve 
% (ar komandu linsolve iegūto rezultātu apzīmēsim ar Xs). 
% Aprēķināt ||X(29) ||2- ||Xs||2.
%  Atbildi dod ar četriem cipariem aiz komata.

A = [8.2,1,4.3,8.6;
     1,10,14.1,2.2;
     4.3,14.1,90,4;
     8.6,2.2,4,31];
B = [-5.5;-3.3;-7.7;9.9];
x = zeros(4,1);

% Оптимальный параметр tau
lambda = eig(A);
tau_opt = 2 / (max(lambda) + min(lambda));

% 29 итераций простой итерации
for k = 1:29
    x = x - tau_opt * (A*x - B);
end

x29 = x;
Xs = linsolve(A,B);

% Разница норм
result = norm(x29,2) - norm(Xs,2);
fprintf('||X(29)|| - ||Xs|| = %.4f\n', result);

%
%Pielietot vienkāršo iterāciju metodi un izpildīt 8 iterācijas, 
%izmantojot sākuma tuvinājumu X = (0000)^T. Izmantot iterāciju parametru t=0.015..
%Kāda ir nesaistes norma ||B – AX(8))||2 ,8. iterācijā? 
%Atbildi dod ar četriem cipariem aiz komata. Вот мой код: 
A = [3,1,9,6;1,10,7,8;9,7,100,13;6,8,13,20;]; 
B = [-3;-2.5;4.7;-3.3;]; 
x = [0;0;0;0]; 
t = 0.015 
for k = 1:8 
    x = x-t*(A*x-B);%vienk.iter.met 
end 
x8 = x; 
result = norm(B-A*x8,2); 
fprintf('%.4f\n',result);

%Dota vienādojumu sistēma AX=B,
%  kur A = [5.2,1.5,9,6.3;1.5,10,10,8.2;9,10,100,13.1;6.3,8.2,13.1,22;];
%  B = [3.5;-3.2;4.5;10]; x = [x1;x2;x3;x4];
%  Pieņemot, ka sistēmu risina ar vienkāršo iterāciju metodi, 
% noskaidrot, par kādu skaitli ir jābūt iterāciju parametra tau vērtībai, 
% lai metode konvergence būtu visātrākā? 
% Atbildi dod ar četriem cipariem aiz komata.

A = [5.2,1.5,9,6.3;
     1.5,10,10,8.2;
     9,10,100,13.1;
     6.3,8.2,13.1,22];
B = [3.5; -3.2; 4.5; 10];

% Находим собственные значения матрицы A
lambda = eig(A);

% Оптимальный параметр τ
tau_opt = 2 / (max(lambda) + min(lambda));

fprintf('Optimālais τ = %.4f\n', tau_opt);
