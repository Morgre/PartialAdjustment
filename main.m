% Author: Moritz Greiwe
% Description: Main File to solve the model on partial adjustment described
% in Bongini (2018)

%% Housekeeping
clear all
clc

%% Parameter initalization
global chi delta r xi tau_d tau_c A s beta

chi = 0.25;
delta = 1;
r = 0.1;
xi = 0.2;
tau_d = 0.2;
tau_c = 0.2;
A = 1;
s = 0.5;
beta = 0.9
%% create grid of shocks + Transition Matrix
z = transpose(linspace(0.5, 1.5, 3));
%P = load("TransMat.mat",'-mat'); for length(z)=40, check for smaller grid
%first
%TransMat = P.P;
P = [0.9, 0.05, 0.05; 0.05, 0.9, 0.05; 0.05, 0.05, 0.9];

%% create grid for capital
k_ss = (chi*A+z(1)/delta)^(1-chi);
k_grid = linspace(0.5*k_ss, 1.5*k_ss, 7);
knext_grid = linspace(0.5*k_ss, 1.5*k_ss, 7);

%% create grid for debt
p_ss = s*(1-delta)*k_ss+A*z(1)*k_ss^chi;
p_grid = linspace(0, 2*p_ss,7);
pnext_grid = linspace(0,2*p_ss,7);


%% Combine grids to grid matrix
k_column = kron(transpose(k_grid), ones(length(knext_grid)*length(p_grid)...
           *length(pnext_grid)*length(z),1));
p_column = kron(ones(length(k_grid),1),kron(transpose(p_grid),...
           ones(length(knext_grid)*length(pnext_grid)*length(z),1)));
knext_column = kron(ones(length(p_grid)*length(k_grid),1),... 
               kron(transpose(knext_grid), ones(length(pnext_grid)*length(z),1)));
pnext_column = kron(ones(length(k_grid)*length(knext_grid)*...
               length(p_grid),1), kron(transpose(pnext_grid), ones(length(z),1)));
z_column = kron(ones(length(k_grid)*length(knext_grid)*length(p_grid)...
           *length(pnext_grid),1), z);
grid_matrix = [k_column, p_column, knext_column, pnext_column, z_column];
%% create transition vector of shocks to calculate expectations
%trMa = kron(P, eye(length(k)*length(k)));
%% Solve constraint problem on the grids using Value function iteration
dev = 1;
maxdev = 10^(-5); % Stopping condition
maxiter = 100000;
iter = 0;
V = 50*ones(length(k_grid)*length(p_grid)*length(z),1);
utility = zeros(length(grid_matrix(:,1)),1);
for i = 1:length(grid_matrix(:,1))
    utility(i,1)= util(grid_matrix(i,1), grid_matrix(i,2),...
                  grid_matrix(i,3), grid_matrix(i,4), grid_matrix(i,5));
end
while (max(abs(dev)) > maxdev) && (iter <= maxiter)
    V_help = reshape(V,length(z), length(V)/length(z));
    V_inter = P*V_help;
    V = reshape(V_inter, length(V),1);
    V_util = beta*kron(ones(length(k_grid)*length(p_grid),1),V)+utility;
    V_util = reshape(V_util, length(z), length(k_grid)*length(p_grid), length(knext_grid)*length(pnext_grid));
    V_util = max(V_util,[], 2);
    Vnew = reshape(V_util, length(k_grid)*length(p_grid)*length(z),1);
    dev = Vnew-V;
    iter = iter+1;
    V = Vnew;
end
