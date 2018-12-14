clear all
tic

% Initial guesses
omega_p_guess = 100;
r_p_guess = 100;
r_d_guess = 100;
x0 = [omega_p_guess r_p_guess r_d_guess ];

% Upper and Lower Bounds
lb = [0, 0.05, 0.2805];
up = [250, 0.1, 0.4];

% Non-linear Programming Solver fmincon
options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
xopt = fmincon(@objective,x0,[],[],[],[],lb,up,@constraint,options)

% Optimisation Results
av_d_opt = calc_av_d(xopt)
size_opt = calc_size(xopt)

% Time to solve
toc 
t = toc

% Define function to calculate velocity of belt
function av_d = calc_av_d(x)
    av_p = x(1);
    r_p = x(2);
    r_d = x(3);
    av_d = (av_p*r_p)/r_d;
end


% Define function to calculate length
function size = calc_size(x)
    av_p = x(1);
    r_p = x(2);
    r_d = x(3);
    size = sqrt(r_p^2+r_p^2) + sqrt(r_d^2+r_d^2) + r_p + r_d;
end


% Define objective function
function obj = objective (x)
    obj = -calc_av_d(x);
end


% Define constraint for optimization
function [c, ceq] = constraint(x)
    c = calc_size(x) - 0.85;
    ceq = [];
end
 