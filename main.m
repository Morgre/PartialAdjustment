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
z = linspace(0.5, 1.5, 40);
P = load("TransMat.mat",'-mat');
TransMat = P.P;
clear P;

%% Solve constraint problem on the grid

dev = 1;
while abs(dev) > maxdev
    Vnew = 
