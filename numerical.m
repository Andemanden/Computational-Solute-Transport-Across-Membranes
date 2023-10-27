clear
close all
clc

% Define parameters
domain_length = 0.5;   % Domain length (meters)
domain_steps = 100;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1);

time_length = 0.5;      % Time length (seconds)
time_steps = 500;       % Number of time steps
dt = time_length / time_steps;

D = 0.01;               % Diffusion coefficient
velocity = 0.3;         % Constant velocity (m/s)
feed_conc = 2; % Constant solute concentration at the first cell
rejection_rate = 0.5; % No rejection

%Constants
R=8.31415; % Gasconstant []
T=273.15+20; %Temperature [K]

% Create a grid for space and time
x = linspace(0, domain_length, domain_steps);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array (1D)
C = ones(domain_steps, time_steps);
% Initial condition (as a column vector)

rho = velocity*dx/dt;
r = D*dx/dt;
fprintf(' r = %f \n rho = %f ', r, rho);

% Time-stepping loop
for j = 2:time_steps
    for i = 1:domain_steps
        % Calculate the second derivative in x direction
        if i == 1
             C(1,:) = feed_conc; % Fill the first cell with inflow concentration
            % No left neighbor at the first cell
            d2Cdx2 = D * (C(2, j-1) -  C(1, j-1)) / dx^2;
             % Apply advection and inflow at the first cell
            dCdt_advection = -velocity * C(i, j-1) / dx;
        elseif i == domain_steps
            % No right neighbor at the last cell and MEMBRANE
            d2Cdx2 = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1)) / dx^2) + C(i,j-1)*(rejection_rate)*velocity/(dx);
        else
            % Calculate the second derivative normally
            d2Cdx2 = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2;
            % Apply advection term
            dCdt_advection = -velocity * (C(i, j-1) - C(i-1, j-1)) / dx;
        end
        % Apply the diffusion-advection-(electromigration) equation
        C(i, j) = C(i, j-1) + dt * (d2Cdx2 + dCdt_advection);
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

% Set axis limits to start at the origin
xlim([0, domain_length]);
ylim([0, time_length]);
zlim([0, max(C(:))]); % Assuming max(C(:)) is the maximum concentration in your data


set(h,'LineStyle','none')