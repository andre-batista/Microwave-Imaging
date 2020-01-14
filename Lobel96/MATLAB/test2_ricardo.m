clc
clear all
load('g.mat');
load('c.mat');

D = C;%(1:2304,1:2304);
I = eye(size(D));
G = gd;

%% 00 inverse algorithm
tic
a = inv(I-G*D);
toc

%% 01 second order approximation - andre18
tic
GD = G*D;
b = I+GD+(GD)^2;
toc

%% 02 c-ricardo
gg = G*G;
tic
GD = G*D;
c = I+GD+D*gg*D;
toc

%% 03 d-ricardo
dig = diag(diag(G));

gg = G*dig;
tic
GD = G*D;
D2 = D*D;
d = I+GD+gg*D2;
toc

%% 04 e-ricardo
G_dig = G-dig; % G - diagonal(G)
didi = dig*dig; % produto das diagonais de G
diGdi = dig*G_dig; % diagonal * (G - diagonal)
Gdidi = G_dig*dig;
Gdi2 = G_dig.^2; % (G - diagonal).*(G - diagonal)

tic
GD = G*D;
D2 = D*D;
%D3 = D2*D;
aux = D*diGdi*D + didi*D2 + Gdidi*D2;
e = I+GD+aux;
toc

%% 05 f-ricardo
tic
GD = G*D;
D2 = D*D;
diD = diag(D);
aux =  diag((Gdi2*diD).*diD)  + D*diGdi*D + didi*D2 + Gdidi*D2;
f = I+GD+aux;
toc

%% tests
norma = norm(a);
%norm(b)
%norm(c)
%(a-di)*d*(a-di)*d + di*d*(a-di)*d + di*d*di*d + (a-di)*d*di*d
%                    d*di*(a-di)*d + di*di*d*d + (a-di)*di*d*d
norm(b-a)/norma
norm(c-a)/norma
norm(d-a)/norma
norm(e-a)/norma
norm(f-a)/norma

fprintf('norma diag')

norma = norm(diag(a));

norm(diag(b-a))/norma
norm(diag(c-a))/norma
norm(diag(d-a))/norma
norm(diag(e-a))/norma
norm(diag(f-a))/norma