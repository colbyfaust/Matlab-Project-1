% Project #1
% Colby Faust, Fisher Barnard, Cameron Mincin, Adam Sabbaghian,
% Louis Veillion
% ME 2543 - Simulations Methods
% Spring 2023

%% Question 3
clc;
close all;
clear;

% lengths given in the tables
L1 = 20;
L2 = 10;
L3 = 30;

% diameter given in the tables
D1 = .167;
D2 = .2083;
D3 = .125;

% Given values of ro, mu, and e from table
ro = .0631;
mu = 0.0000205;
e = .00085;

% converting the diameter to area
A1 = pi*(D1/2)^2;
A2 = pi*(D2/2)^2;
A3 = pi*(D3/2)^2;

% range of flow rates
Q_range = 100:50:1500;

% arrays to store computed values of Q and deltaP
Q_array = zeros(size(Q_range));
deltaP_array = zeros(size(Q_range));

% options for fsolve
options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',1e4,'MaxIterations',1e4);

for i = 1:length(Q_range)
    
    % initial guess for Q1, Q2, Q3, and deltaP
    Q0 = [160; 176; 137];
    deltaP0 = 4;
    
    % total flow rate Q in gpm
    Q = Q_range(i);
    
    % function for fsolve
    fun = @(x) [x(1)+x(2)+x(3)-Q;
                64/((ro*x(1)/A1)*D1/mu)*(L1/D1)*((x(1)/A1)^2/2)-x(4)/ro;
                0.25*(log(((e/D2)/3.7)+5.74/((ro*x(2)/A2)*D2/mu))^(-2))*(L2/D2)*((x(2)/A2)^2/2)-x(4)/ro;
                0.25*(log(((e/D3)/3.7)+5.74/((ro*x(3)/A3)*D3/mu))^(-2))*(L3/D3)*((x(3)/A3)^2/2)-x(4)/ro];

    % solve for Q1, Q2, Q3, and deltaP
    x = fsolve(fun, [Q0; deltaP0], options);

    % store the computed values of Q and deltaP
    Q_array(i) = Q;
    deltaP_array(i) = x(4);
end

% plot the variation of Q with deltaP
plot(deltaP_array, Q_array);
xlabel('Change in pressure (psi)');
ylabel('Total flow rate Q (gpm)');
title('Total flow rate Vs. Change in pressure');
