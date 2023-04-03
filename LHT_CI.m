% Modify this as needed :---------------
alfa = 0.05 %Change parameter as needed
fnam = "test_data.csv"; % Change file name
% Change directory path here, live commented if file is in current directory:
% dir_path = "/Users/datasets/"

% Begin calculations :---------------
Z = norminv(1-alfa/2);
T = readtable(fnam)
t = T.t;
n = T.n;
h = T.h;
K = sum(h);
N = n(1);
g = h/K;
x = n(1:end-1)-n(2:end);
x=[x;0];
f = x/sum(x);

% Table to use:

new_table =[t n h g f];
T = table(t, n, h, g, f, 'VariableNames', {'t', 'n', 'h', 'g', 'f'});
fnam2 = strrep(fnam, '.csv', '_added.csv');
writetable(T, fnam2);

% 1) Longevity
EX = t'*f;
EX2=(t.*t)'*f;
VX = EX2-(EX*EX);
%[EX VX]

% CI for longevity:
CI_L = [EX VX EX-Z*sqrt(VX/N) EX+Z*sqrt(VX/N)];

% 2) R0

R0 = K/N; %(from paper)

% 3) Generation time

ET = t'*g;
ET2=(t.*t)'*g;
VT = ET2-(ET*ET);
%[ET VT]

% CI for longevity:
CI_G = [ET VT ET-Z*sqrt(VT/K) ET+Z*sqrt(VT/K)];

% 4) Growth rate r

%Define function handle
f = @(x) lambd(x,t,g,1/R0);
% initial value is lambdac
r0 = 0.02;
r = fzero(f,r0);

mu =1/R0;
s2 = exp(-2*r*t')*g-mu*mu;
c = Z*sqrt(s2/N);

f = @(x) lambd(x,t,g,mu+c);
r_low = fzero(f,r);

f = @(x) lambd(x,t,g,mu-c);
r_upp = fzero(f,r);
rr = [r r_low r_upp];

lam = exp(rr);
clc
disp(sprintf('Initial number of individuals N : %s', num2str(N)))
disp(sprintf('Offspring size K                : %s', num2str(K)))
disp(sprintf('R0                              : %s', num2str(R0)))
disp(sprintf('Longevity                       : %s', num2str(CI_L)))
disp(sprintf('Generation time                 : %s', num2str(CI_G)))
disp(sprintf('r                               : %s', num2str(rr)))
disp(sprintf('lambda                          : %s', num2str(lam)))

disp(sprintf('New data saved to: %s', fnam2))



function out = lambd(r,t,g,W)
out = sum(exp(-r*t).*g)-W;
end



