% MODIFY PARAMETERS AS NEEDED
alfa = 0.05; % Change parameter as needed
fnam = "test_data_18.csv"; % Change file name

% IMPORT DATA
T = readtable(fnam);
t = T.t;
n = T.n;
h = T.h;

% CALCULATE VALUES
Z = norminv(1-alfa/2);
K = sum(h);
N = n(1);
g = h/K;
x = n(1:end-1) - n(2:end);
x = [x; 0];
f = x/sum(x);

% CREATE A NEW TABLE
new_table = [t, n, h, g, f];
T = table(t, n, h, g, f, 'VariableNames', {'t', 'n', 'h', 'g', 'f'});
fnam2 = strrep(fnam, '.csv', '_added.csv');
writetable(T, fnam2);

% CALCULATE LONGEVITY
EX = t' * f;
EX2 = (t .* t)' * f;
VX = EX2 - (EX * EX);
CI_L = [EX, VX, EX - Z * sqrt(VX/N), EX + Z * sqrt(VX/N)];

% CALCULATE R0
R0 = K / N;

% CALCULATE GENERATION TIME
ET = t' * g;
ET2 = (t .* t)' * g;
VT = ET2 - (ET * ET);
CI_G = [ET, VT, ET - Z * sqrt(VT/K), ET + Z * sqrt(VT/K)];

% CALCULATE GROWTH RATE r
f = @(x) lambd(x, t, g, 1/R0);
r0 = 0.0;
r = fzero(f, r0);

mu = 1/R0;
s2 = exp(-2 * r * t') * g - mu * mu;
c = Z * sqrt(s2/N);

f = @(x) lambd(x, t, g, mu + c);
r_low = fzero(f, r);

f = @(x) lambd(x, t, g, mu - c);
r_upp = fzero(f, r);
rr = [r, r_low, r_upp];

lam = exp(rr);
clc

% DISPLAY RESULTS
disp(sprintf('Initial number of individuals N : %s', num2str(N)))
disp(sprintf('Offspring size K                : %s', num2str(K)))
disp(sprintf('R0                              : %s', num2str(R0)))
disp(sprintf('Longevity                       : %s', num2str(CI_L)))
disp(sprintf('Generation time                 : %s', num2str(CI_G)))
disp(sprintf('r                               : %s', num2str(rr)))
disp(sprintf('lambda                          : %s', num2str(lam)))
disp(sprintf('New data saved to: %s', fnam2))

% LAMBDA FUNCTION
function out = lambd(r, t, g, W)
    out = sum(exp(-r * t) .* g) - W;
end
