clear
close all
clc

% Define parameters
domain_length = 0.5;   % Domain length (meters)
discretization = 50;    % Number of spatial discretization points
dx = domain_length / (discretization - 1);

time_length = 5;      % Time length (seconds)
time_steps = 5000;       % Number of time steps
dt = time_length / time_steps;

D = 0.001;               % Diffusion coefficient
velocity = 0.1;         % Constant velocity (m/s)
inflow_concentration = 2.0; % Constant solute concentration at the first cell


% Create a grid for space and time
x = linspace(0, domain_length, discretization);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array (1D)
C = ones(discretization, time_steps);

% Set up the diffusion-convection ODE
rejection_rate = 0; % No rejection

% Initial condition (as a column vector)
C(:, 1) = zeros(discretization, 1); % Initialize to 0 everywhere
C(1, 1) = inflow_concentration; % Set the inflow concentration at the first position

% Time-stepping loop
for k = 2:time_steps
    for i = 1:discretization
        % Calculate the second derivative in x direction
        if i == 1
            % No left neighbor at the first cell
            d2Cdx2 = D * (C(2, k-1) -  C(1, k-1)) / dx^2;
             % Apply advection and inflow at the first cell
            dCdt_convection = -velocity * C(i, k-1) / dx;
            C(1, k) = inflow_concentration; % Fill the first cell with inflow concentration
        elseif i == discretization
            % No right neighbor at the last cell
            d2Cdx2 = D * (C(discretization - 1, k-1) - C(discretization, k-1)) / dx^2;
        else
            % Calculate the second derivative normally
            d2Cdx2 = D * (C(i + 1, k-1) - 2 * C(i, k-1) + C(i - 1, k-1)) / dx^2;
            % Apply advection term
            dCdt_convection = -velocity * (C(i, k-1) - C(i-1, k-1)) / dx;
        end
        % Apply the diffusion-convection equation
        C(i, k) = C(i, k-1) + dt * (d2Cdx2 - rejection_rate * C(i, k-1) + dCdt_convection);
    end
end


% Define fractions of time steps you want to visualize
time_fraction = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];  % For example, 0.1 corresponds to 10% of time steps

% Calculate the corresponding time indices
time_instances = round(time_fraction * time_steps);

% Create a figure for the first plot
figure;
hold on;

% Plot concentration at specific time instances
for i = 1:length(time_instances)
    time_index = time_instances(i);
    plot(x, C(:, time_index));
end

xlabel('Position (meters)');
ylabel('Concentration');
title('Concentration Over Time at Different Instances');

% Add a legend for clarity
legend(arrayfun(@(f) ['t=', num2str(f)], time_fraction, 'UniformOutput', false));

grid on;
hold off;

% Create a 3D surface plot to visualize concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, C); % Transpose removed here
xlabel('Position (meters)');
ylabel('Time (seconds)');
zlabel('Concentration');
title('Concentration Over Time and Position');

set(h,'LineStyle','none')