% Simple example finding all solution to a synthetic problem
clear all; close all;
addpath('..');

%% Generate data
rows = 5;
cols = 5;
num_nodes = rows*cols;
disc = 10^2;

% E = unary + \lambda_1*slope1 + \lambda_2*slope2 + pairwise
slope1 = round((randn(num_nodes,1))*disc);
slope2 = round((randn(num_nodes,1))*disc);
lambda_0 = round(randn(size(slope1))*disc);
unary = - lambda_0.*slope1;

% Generate pairwise weights
% 4 connectivity neighborhood
delta_connectivity = [-1 0; 0 -1; 1 0; 0 1]; 

% Generate all pairs
connectivity = Parametric.create_connectivty(rows, cols, delta_connectivity);

% Set weights
pairwise = [connectivity round((abs(randn(size(connectivity,1),1)))*disc) zeros(size(connectivity,1),1)];

% Create object 
P = Parametric(unary, slope1, slope2, pairwise);
P.iterations = 5000;

%%
% Evaluate one point (\lambda, \mu)
[one_solution, energy] = P.evaluate(0,0);

% Solve
solutions = P.find_all_solutions();
Parametric.draw_solution_diagram(solutions);