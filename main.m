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
beta = 0.9;
%% create grid of shocks + Transition Matrix
z = transpose(linspace(0.5, 1.5, 3));
%P = load("TransMat.mat",'-mat'); for length(z)=40, check for smaller grid
%first
%TransMat = P.P;
P = [0.5, 0.25, 0.25; 0.25, 0.5, 0.25; 0.25, 0.25, 0.5];

%% create grid for capital
k_ss = (chi*A+z(1)/delta)^(1-chi);
k_grid = linspace(0.5*k_ss, 1.5*k_ss, 31);
knext_grid = linspace(0.5*k_ss, 1.5*k_ss, 31);

%% create grid for debt
p_ss = s*(1-delta)*k_ss+A*z(1)*k_ss^chi;
p_grid = linspace(-2*p_ss, 2*p_ss,31);
pnext_grid = linspace(-2*p_ss,2*p_ss,31);


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

%% Solve constraint problem on the grids using Value function iteration
dev = 1;
maxdev = 10^(-5); % Stopping condition
maxiter = 10000;
iter = 0;
V = ones(length(k_grid)*length(p_grid)*length(z),1);
utility = util(grid_matrix(:,1),grid_matrix(:,2),grid_matrix(:,3),grid_matrix(:,4),grid_matrix(:,5));
while (max(abs(dev)) > maxdev) && (iter <= maxiter)
    V_help = reshape(V,length(z), length(V)/length(z));
    V_inter = kron(ones(length(k_grid)*length(p_grid),1),P*V_help);
    V_inter = reshape(V_inter, numel(V_inter),1); 
    V_util = beta*V_inter+utility;
    V_util = reshape(V_util, length(z), length(k_grid)*length(p_grid), length(knext_grid)*length(pnext_grid));
    [V_util2, index] = max(V_util,[], 2);
    Vnew = reshape(V_util2, length(k_grid)*length(p_grid)*length(z),1);
    dev = Vnew-V;
    iter = iter+1;
    V = Vnew;
end

%% Plot policy functions
index = reshape(index, length(k_grid)*length(p_grid)*length(z),1);
policy_plot_z1 = reshape(index(1:3:end), length(p_grid), length(k_grid)); 
policy_plot_z2 = reshape(index(2:3:end), length(p_grid), length(k_grid));
policy_plot_z3 = reshape(index(3:3:end), length(p_grid), length(k_grid));

subplot(1,3,1);
surf(policy_plot_z1);
title('Policy function for debt, Shock low')
xlabel('Debt grid point')
ylabel('Capital grid point')
zlabel('Debt tomorrow')

subplot(1,3,2);
surf(policy_plot_z2);
title('Policy function for debt, Shock mid')
xlabel('Debt grid point')
ylabel('Capital grid point')
zlabel('Debt tomorrow')

subplot(1,3,3);
surf(policy_plot_z3);
title('Policy function for debt, Shock high')
xlabel('Debt grid point')
ylabel('Capital grid point')
zlabel('Debt tomorrow')