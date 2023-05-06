% Project #1
% Colby Faust, Fisher Barnard, Cameron Mincin, Adam Sabbaghian,
% Louis Veillion
% ME 2543 - Simulations Methods
% Spring 2023

%% Question 1
clc;
close all;
clear;

% lengths given in the tables
L1 = 10;
L2 = 10;
L3 = 10;

% diameter given in the tables
D1 = .167;
D2 = .167;
D3 = .167;

% Given values of ro, mu, and e from table
ro = .0631;
mu = 0.0000205;
e = .00085;

% converting the diameter to area
A1 = pi*(D1/2)^2;
A2 = pi*(D2/2)^2;
A3 = pi*(D3/2)^2;

% initial guess for Q1, Q2, Q3, and deltaP
Q0 = [120; 150; 200];
deltaP0 = 4;

% total flow rate Q = 750 gpm in L/s
Qf_1 = 47.3176;

% function for fsolve
fun = @(x) [x(1)+x(2)+x(3)-47.3176;
            0.25*(log(((e/D1)/3.7)+5.74/((ro*x(2)/A2)*D1/mu))^(-2))*(L1/D1)*((x(1)/A1)^2/2)-x(4)/ro;
            0.25*(log(((e/D2)/3.7)+5.74/((ro*x(2)/A2)*D2/mu))^(-2))*(L2/D2)*((x(2)/A2)^2/2)-x(4)/ro;
            0.25*(log(((e/D3)/3.7)+5.74/((ro*x(3)/A3)*D3/mu))^(-2))*(L3/D3)*((x(3)/A3)^2/2)-x(4)/ro];

% options for fsolve
options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e4,'MaxIterations',1e4);

% fsolve to give us values for Q1 Q2 Q3 and delta P
x = fsolve(fun, [Q0; deltaP0], options);

% converting the values of Q's from L/s to gpm
Q1 = (x(1))*15.85032;
Q2 = (x(2))*15.85032;
Q3 = (x(3))*15.85032;
Qf_2 = 47.3176*15.85032;

% converting deltaP from N/m^2 to lbf/ft^2
deltaP = (x(4))*.0208854;


% fprintf function to display the values of the unknowns
fprintf('The value of Q1 is %.4f gpm\n',Q1);
fprintf('The value of Q2 is %.4f gpm\n',Q2);
fprintf('The value of Q3 is %.4f gpm\n',Q3);
fprintf('The value of Q is %.4f gpm\n',Qf_2);
fprintf('The value of deltaP is %.4f lbf/f^2\n',deltaP);

