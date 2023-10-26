clear all
close all
clc
tic
 
%hej
% Define parameters
domain_length = 0.5;       % Domain length (meters)
discretization = 50;       % Number of spatial discretization points
dx = domain_length / (discretization - 1);

time_length = 0.5;         % Time length (seconds)
time_steps = 500;          % Number of time steps
dt = time_length / time_steps;

D = 0.01;                   % Diffusion coefficient

% Create a grid for space and time
x = linspace(0, domain_length, discretization);
t = linspace(0, time_length, time_steps);

% Set up the diffusion ODE
diffusion_ode = @(t, C) D * second_derivative(C, dx);

% Initial condition (as a column vector)
C0 = ones(discretization, 1); % Initialize to 1 everywhere
C0(1) = 2; % Set the concentration at the first position to 2

% Solve the diffusion ODE
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, C] = ode15s(diffusion_ode, t, C0, options);

% Extract concentration at a specific time
time_index = 10;  % Change this index to select the time you want to visualize
concentration_at_specific_time = C(time_index, :);

% Plot concentration over distance at a specific time
figure;
plot(x, concentration_at_specific_time);
xlabel('Position (meters)');
ylabel('Concentration');
title(['Concentration Over Distance at Time: ' num2str(t(time_index))]);

% Create a 3D surface plot to visualize concentration over time
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, C'); % Transpose C for correct dimensions
xlabel('Position (meters)');
ylabel('Time (seconds)');
zlabel('Concentration');
title('Concentration Over Time and Position');
% Set axis limits to start at the origin
xlim([0, domain_length]);
ylim([0, time_length]);
zlim([0, max(C(:))]); % Assuming max(C(:)) is the maximum concentration in your data

set(h,'LineStyle','none')
% Function to compute the second derivative using finite differences
function d2Cdx2 = second_derivative(C, dx)
    d2Cdx2 = (C([2:end, end]) - 2 * C + [C(1); C(1:end-1)]) / dx^2; % Make it a column vector
end
