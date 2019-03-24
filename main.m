% Author: Moritz Greiwe
% Description: Main File to solve the model on partial adjustment described
% in Bongini (2018)
global chi delta r xi tau_d tau_c A s

%% Parameter initalization
chi = 0.25;
delta = 1;
r = 1;
xi = 0.2;
tau_d = 0.2;
tau_c = 0.2;
A = 1;
s = 0.5;

%% create grid of shocks + Transition Matrix
z = transpose(linspace(0.5, 1.5, 5));
%P = load("TransMat.mat",'-mat'); for length(z)=40, check for smaller grid
%first
%TransMat = P.P;
P = [0.5, 0.25, 0.25; 0.25, 0.5, 0.25; 0.25, 0.25, 0.5];

%% create grid for capital
k_ss = 1;
k_grid = linspace(0.5*k_ss, 1.5*k_ss, 7);
knext_grid = linspace(0.5*k_ss, 1.5*k_ss, 7);

%% create grid for debt
p_ss = s*(1-delta)*k_ss+A*z(1)*k_ss^chi;
p_grid = linspace(0, 2*p_ss,7);
pnext_grid = linspace(0,2*p_ss,7);


%% Combine grids to grid matrix
k_column = kron(transpose(k_grid), ones(length(knext_grid)*length(p_grid)*length(pnext_grid)*length(z),1));
p_column = kron(kron(transpose(p_grid), ones(length(knext_grid)*length(pnext_grid)*length(z),1)),ones(length(k_grid),1));
knext_column = kron(kron(transpose(knext_grid), ones(length(pnext_grid)*length(z),1)), ones(length(p_grid)*length(k_grid),1));
pnext_column = kron(kron(transpose(pnext_grid), ones(length(z),1)),ones(length(k_grid)*length(knext_grid)*length(p_grid),1));
z_column = kron(ones(length(k_grid)*length(knext_grid)*length(p_grid)*length(pnext_grid),1), z);
grid_matrix = [k_column, p_column, knext_column, pnext_column, z_column];
%% create transition vector for states
%trMa = kron(P, eye(length(k)*length(k)));
%% Solve constraint problem on the grids using Value function iteration
dev = 1;
maxdev = 10^(-5); % Stopping condition
maxiter = 10000;
iter = 0;
V = k_ss*(length(k)^2*length(z));
while (abs(dev) > maxdev) && (iter <= maxiter)
    utility = zeros(length(grid_matrix(:,1)),1);
    for i = 1:length(grid_matrix(:,1))
        utility(i,1)= util(grid_matrix(i,1), grid_matrix(i,2), grid_matrix(i,3), grid_matrix(i,4), grid_matrix(i,5));
    end
%    Vnew = utility + trMa*V; %include transition matrix between the states                                                                                                          
%    dev = Vnew-V;
%    iter = iter+1;
%    V = Vnew;
end
