clear
close all
clc

% Define parameters
domain_length = 0.5;   % Domain length (meters)
discretization = 50;    % Number of spatial discretization points
dx = domain_length / (discretization - 1);

time_length = 0.5;      % Time length (seconds)
time_steps = 500;       % Number of time steps
dt = time_length / time_steps;

D = 0.01;               % Diffusion coefficient
velocity = 0.1;         % Constant velocity (m/s)
inflow_concentration = 1.0; % Constant solute concentration at the first cell

% Create a grid for space and time
x = linspace(0, domain_length, discretization);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array (1D)
C = zeros(discretization, time_steps);

% Set up the diffusion-convection ODE
rejection_rate = 0; % No rejection

% Initial condition (as a column vector)
C(:, 1) = zeros(discretization, 1); % Initialize to 0 everywhere
C(discretization/2, 1) = inflow_concentration; % Set the inflow concentration at the first position

% Time-stepping loop
for k = 2:time_steps
    for i = 1:discretization
        % Calculate the second derivative in x direction
        d2Cdx2 = D * (C(min(i + 1, discretization), k-1) - 2 * C(i, k-1) + C(max(i - 1, 1), k-1)) / dx^2;
        
        % Calculate convection term
        dCdt_convection = -velocity * (C(i, k-1) - C(max(i - 1, 1), k-1)) / dx;
        
        % Apply the diffusion-convection equation
        C(i, k) = C(i, k-1) + dt * (d2Cdx2 - rejection_rate * C(i, k-1) + dCdt_convection);
    end
end

% Plot concentration at a specific time
time_index = 10;  % Change this index to select the time you want to visualize
figure;
plot(x, C(:, time_index));
xlabel('Position (meters)');
ylabel('Concentration');
title(['Concentration Over 1D Domain at Time: ' num2str(t(time_index))]);


% Create a 3D surface plot to visualize concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, C); % Transpose removed here
xlabel('Position (meters)');
ylabel('Time (seconds)');
zlabel('Concentration');
title('Concentration Over Time and Position');

set(h,'LineStyle','none')