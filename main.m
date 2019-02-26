% Author: Moritz Greiwe
% Description: Main File to solve the model on partial adjustment described
% in Bongini (2018)
global chi delta r xi tau_d tau_c A

%% Parameter initalization
chi = 0.25;
delta = 1;
r = 1;
xi = 0.2;
tau_d = 0.2;
tau_c = 0.2;
A = 1;

%% create grid of shocks + Transition Matrix
z = transpose(linspace(0.5, 1.5, 3));
%P = load("TransMat.mat",'-mat'); for length(z)=40, check for smaller grid
%first
%TransMat = P.P;
P = [0.5, 0.25, 0.25; 0.25, 0.5, 0.25; 0.25, 0.25, 0.5];

%% create grid for capital
k_ss = 1
k_grid = linspace(0.5*k_ss, 1.5*k_ss, 5);

%% Combine grids to create grid-vectors
%k_grid_vector = zeros(length(k_grid)*length(z),1);
%for i=1:length(k_grid)
%    down = (i-1)*length(z)+1;
%    up = i*length;
%    k_grid_vector(down:up,1) = ones(length(z))*k_grid(i);
%end
%z_grid = kron(ones(length(k_grid),1), z);

%% create transition vector for states
trMa = kron(P, eye(length(k)*length(k)));
%% Solve constraint problem on the grids using Value function iteration
dev = 1;
maxdev = 10^(-5); % Stopping condition
maxiter = 10000;
iter = 0;
V = k_ss*(length(k)^2*length(z));
while (abs(dev) > maxdev) && (iter <= maxiter)
    Vnew = util(z, k_grid, k_grid) + trMa*V; %include transition matrix between the states                                                                                                          
    dev = Vnew-V;
    iter = iter+1;
    V = Vnew;
end
