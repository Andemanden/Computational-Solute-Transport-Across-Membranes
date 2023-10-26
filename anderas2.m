clear
close all
clc

% Define parameters
domain_length = 0.5;       % Domain length (meters)
discretization = 50;       % Number of spatial discretization points
dx = domain_length / (discretization - 1);

time_length = 0.5;         % Time length (seconds)
time_steps = 500;          % Number of time steps
dt = time_length / time_steps;

D = 0.01;                   % Diffusion coefficient
membrane_position = domain_length; % Position of the membrane (meters)
rejection_rate = 1;     % Rejection rate (e.g., 20% rejection)

% Create a grid for space and time
x = linspace(0, domain_length, discretization);
t = linspace(0, time_length, time_steps);

% Initialize the concentration matrix
C = zeros(discretization, time_steps);

% Set up the diffusion ODE
diffusion_ode = @(t, y) D * second_derivative(y, dx, D);

% Initial condition (as a column vector)
C(:, 1) = ones(discretization, 1); % Initialize to 1 everywhere
C(1, 1) = 2; % Set the concentration at the first position to 2

% Solve the diffusion ODE
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Initialize the rejected material
rejected = zeros(1, time_steps);

% Time-stepping loop
for i = 2:time_steps
    tspan = [t(i-1), t(i)];
    
    % Set the initial condition for the ODE solver
    initial_condition = C(:, i-1);
    
    % Use the ODE solver to solve the diffusion equation
    [~, result] = ode15s(diffusion_ode, tspan, initial_condition, options);
    
    % Reshape the result to match the size of C
    C(:, i) = result(end, :)';
    
    % Calculate the amount rejected at the membrane
    if x(1) < membrane_position
        amount_rejected = rejection_rate * (C(1, i-1) - C(2, i-1));
        % Adjust the concentration to account for rejection at the cell before the membrane
        C(1, i) = C(1, i) - amount_rejected;
    else
        amount_rejected = 0;
    end
    
    % Accumulate rejected material
    rejected(i) = rejected(i-1) + amount_rejected;
end

% Plot concentration over distance at a specific time
time_index = 10;  % Change this index to select the time you want to visualize
figure;
plot(x, C(:, time_index));
xlabel('Position (meters)');
ylabel('Concentration');
title(['Concentration Over Distance at Time: ' num2str(t(time_index))]);

% Create a 3D surface plot to visualize concentration over time
[T, X] = meshgrid(t, x);
figure;
h = surf(T, X, C); % Swap T and X for the surf plot
xlabel('Time (seconds)');
ylabel('Position (meters)');
zlabel('Concentration');
title('Concentration Over Time and Position');

% Set axis limits to start at the origin
xlim([0, time_length]);
ylim([0, domain_length]);
zlim([0, max(C(:))]); % Assuming max(C(:)) is the maximum concentration in your data

set(h, 'LineStyle', 'none');

% Function to compute the second derivative using finite differences
function dydt = second_derivative(y, dx, D)
    % Calculate the second derivative
    dydt = D * (circshift(y, 1) - 2 * y + circshift(y, -1)) / dx^2;
end
